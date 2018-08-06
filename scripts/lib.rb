require_relative('utils')
require('git')
require('fileutils')
require('open3')

SLEEP_TIME = 1 # seconds
MAX_ITER = 4

def is_pull_request?
  travis_pr = ENV['TRAVIS_PULL_REQUEST']
  if travis_pr == "false"
    false
  elsif travis_pr.nil? or travis_pr.empty?
    # OK, we're not on Travis...
    av_pr = ENV['APPVEYOR_PULL_REQUEST_NUMBER']
    if av_pr.nil? or av_pr.empty?
      false
    else
      true
    end
  else
    # TRAVIS_PULL_REQUEST exists and is not false so we *are* on a pull request
    # on Travis; return true
    true
  end
end

# String String String -> Git::Base
# Robustly clones a source repo to the given target
# repo_url: String   = URL to repository to clone from
# local_path: String = File path to where the cloned repository should live on
#                      the local filesystem
# remote_branch: String = string of branch to clone
def robust_clone(repo_url, local_path, remote_branch)
  g = nil
  if File.exist?(local_path)
    if File.exist?(File.join(local_path, '.git'))
      puts("  .git directory exists at local_path")
      g = Git.open(local_path)
      r = g.remotes.select {|r| r.url == repo_url}[0]
      origin = g.remotes.select {|r| r.name == "origin"}[0]
      if r.nil? or r != origin
        puts("  re-specifying origin")
        g.remove_remote("origin") unless origin.nil?
        g.add_remote("origin", repo_url)
      end
      puts("  fetching changes")
      `cd #{local_path} && git fetch`
      if UTILS::git_remote_branch_exists?(remote_url, remote_branch)
        puts("  remote branch (#{remote_branch}) exists, checking out...")
        `cd #{local_path} && git checkout -b #{remote_branch} origin/#{remote_branch}`
      else
        puts("  remote branch doesn't exist, creating local branch with name #{remote_branch}")
        `cd #{local_path} && git checkout -b #{remote_branch}`
      end
    else
      puts("  local_path exists but no .git directory in it...")
      if Dir[File.join(local_path, '*')].empty?
        if UTILS::git_remote_branch_exists?(remote_url, remote_branch)
          puts("  remote branch (#{remote_branch}) exists, cloning it...")
          `cd #{local_path} && git clone -b #{remote_branch} #{repo_url} .`
        else
          puts("  remote branch (#{remote_branch}) does NOT exist, checking out default branch and making new local branch")
          `cd #{local_path} && git clone #{repo_url} . && git checkout -b #{remote_branch}`
        end
        g = Git.open(local_path)
      else
        raise "  local_path \"#{local_path}\" exists but is NOT empty; please remove files and try again"
      end
    end
  else
    FileUtils.mkdir_p(local_path)
    if UTILS::git_remote_branch_exists?(repo_url, remote_branch)
      puts("  Cloning remote_branch: #{remote_branch}")
      cmd = "cd #{local_path} && git clone -b #{remote_branch} #{repo_url} ."
      `#{cmd}`
    else
      puts("  Remote_branch (#{remote_branch}) DOES NOT exist; cloning default and creating local branch")
      `cd #{local_path} && git clone #{repo_url} . && git checkout -b #{remote_branch}`
    end
    g = Git.open(local_path)
  end
  g
end

def find_commit_from_substr(g, substr)
  c = nil
  g.tags.each do |t|
    if t.name.include?(substr)
      c = t.objectish
      break
    end
  end
  c
end

def determine_branch
  if ! ENV['SRC_BRANCH'].nil?
    ENV['SRC_BRANCH']
  elsif ! ENV['TRAVIS_BRANCH'].nil?
    puts("Tried to find environment variable SRC_BRANCH but it was nil...")
    puts("trying TRAVIS_BRANCH...")
    ENV['TRAVIS_BRANCH']
  elsif ! ENV['APPVEYOR_REPO_BRANCH'].nil?
    puts("Tried to find environment variable SRC_BRANCH but it was nil...")
    puts("trying APPVEYOR_REPO_BRANCH...")
    ENV['APPVEYOR_REPO_BRANCH']
  else
    puts("Could not determine branch... exiting with non-zero status")
    exit(1)
  end
end

def mimic_source(g, branch, src_commit)
  puts("  Mimicking source repository")
  puts("  - branch: #{branch}")
  puts("  - src_commit: #{src_commit}")
  g.branch(branch).checkout
end

def run_case(in_root, out_root, arch, a, c)
  path = File.join(out_root, arch, c[:id])
  FileUtils.mkdir_p(path) unless File.exist?(path)
  FileUtils.cp(File.join(in_root, a[:exe]), path)
  c[:input_files].each do |f|
    FileUtils.cp(File.join(in_root, f), path)
  end
  cmd = [
    "cd #{path}",
    "#{a[:cmd]} #{c[:args].join(' ')}"
  ].join(" && ")
  stdout, stderr, exitcode = Open3.capture3(cmd)
  files_to_add = c[:output_files].map {|o|
    File.join(path, File.basename(o))
  }
  w = lambda do |name, data|
    p = File.join(path, name)
    File.write(p, data)
    files_to_add << p
  end
  w['stdout.log', UTILS.universal_newline(stdout)]
  w['stderr.log', UTILS.universal_newline(stderr)]
  #w['exitcode.log', exitcode]
  c[:input_files].each do |f|
    FileUtils.rm(File.join(path, File.basename(f)))
  end
  FileUtils.rm(File.join(path, File.basename(a[:exe])))
  files_to_add
end

# robust pull/push
def robust_push_pull(g, branch, the_commit, the_tag, rt_url, our_dir=nil, to_be_deleted=nil)
  1.upto(MAX_ITER).each do |try_no|
    puts("    Attempt #{try_no} (of #{MAX_ITER})")
    begin
      g.pull(rt_url, branch) if g.is_remote_branch?(branch)
      # -X ignore-space-at-eol
      #cmd = "git pull -X ours --allow-unrelated-histories #{rt_url} #{branch}"
      #cmd = "git pull #{rt_url} #{branch}"
      #puts("running command: #{cmd}")
      #out, err, status = Open3.capture3(cmd) if g.is_remote_branch?(branch)
      #raise 'Error' if status != 0
      g.add(all: true)
      code = system("cd #{g.dir} && git diff --quiet --cached --exit-code")
      if !code
        puts("    Added all")
        g.commit("Committing changes from merge...\n#{the_commit}")
        puts("    Commit successful")
      else
        puts("    No changes to commit moving forward")
      end
    rescue => e
      puts("    Trying to fix suspected auto-merge conflict")
      puts("    Error: #{e.message}")
      #puts("    stdout: #{out}")
      #puts("    stderr: #{err}")
      puts("    Git repo directory: #{g.dir}")
      # OK, we probably have a merge conflict
      # we need to:
      # 1. find which files are in conflict
      # 2. `git checkout --ours #{filename}` for each file in conflict
      # 3. re-add the files
      # 4. re-commit the files with some sort of commit message
      chngs = UTILS.list_changes(g.dir)
      # `git diff --name-only --diff-filter=U`
      puts("    Asking git for list of conflicted files")
      files_in_conflict = `cd #{g.dir} && git diff --name-only --diff-filter=U`
      puts("    files in conflict: #{files_in_conflict.split.inspect}")
      puts("    " + "="*60)
      puts("    diff of conflicts:")
      puts("#{`cd #{g.dir} && git diff --ws-error-highlight=all`}")
      puts("    " + "="*60)
      num_ours = 0
      num_theirs = 0
      files_in_conflict.split.each do |fname|
        puts("checking out and re-adding #{fname}")
        # If under our architecture, use --ours; else --theirs...
        whos = nil
        if our_dir.nil? or fname.include?(our_dir)
          whos = "ours"
        else
          whos = "theirs"
        end
        deleted = false
        d1 = chngs['Deleted'].include?(fname) 
        d2 = (!to_be_deleted.nil? and to_be_deleted.include?(fname))
        if d1 or d2
          `cd #{g.dir} && git rm #{fname}`
          deleted = true
        else
          `cd #{g.dir} && git checkout --#{whos} #{fname}`
        end
        if $?.exitstatus == 0
          g.add(fname) if not deleted
          if whos == "ours"
            num_ours += 1
          else
            num_theirs += 1
          end
        end
      end
      puts("    Number added: theirs: #{num_theirs} ours: #{num_ours}")
      if files_in_conflict.length > 0
        code = system("cd #{g.dir} && git diff --quiet --cached --exit-code")
        if !code
          puts("    (Re-)Committing...")
          g.commit("Commit to fix auto-merge conflict\n#{the_commit}")
        else
          puts("    No changes to commit, moving forward")
        end
      else
        puts("    Number of files in conflict is 0 (warning: how did we get into this branch?)")
      end
      puts("    Done!")
    end
    puts("    Pushing to origin!")
    begin
      g.push(rt_url, "HEAD:#{branch}", {:tags=>true})
    rescue => e
      # Possible that we have an error related to remote tagging
      # Let's check if the word "tag" is in the error message
      puts("    -------")
      puts("    push caused an error...\n#{e.message}")
      if e.message.include?("tag")
        # if tag already exists on remote, we need to delete, force annotate, and push again.
        puts("    apparently, tag already exists on remote!")
        puts("    attempt a retag...")
        retag(g.dir, the_tag, rt_url)
        puts("    back from retag...")
        next
      else
        # We don't know what happened. Report error and bail.
        puts("    Don't know how to handle this error... exiting...")
        puts("    -------")
        exit(1)
      end
    end
    break
  end
  puts("  Done robust_pull_push!")
end

def retag(git_dir, tag_name, rt_url)
  puts("Tag, #{tag_name}, exists")
  1.upto(MAX_ITER).each do
    sleep(SLEEP_TIME)
    # delete tag on remote
    # See https://stackoverflow.com/a/19300065
    puts('    Deleting tag on remote...')
    `cd #{git_dir} && git push --quiet "#{rt_url}" :refs/tags/#{tag_name}`
    puts("    Remote deletion attempted. Return code: #{$?}")
    tag_output = `cd #{git_dir} && git ls-remote "#{rt_url}" refs/tags/#{tag_name}`
    puts("    Tag query attempted. Return value: #{tag_output}")
    tag_deleted = tag_output.strip.empty?
    break if tag_deleted
  end
  # force annotate the tag again
  puts("    Force annotating the local tag again")
  `cd #{git_dir} && git tag -fa #{tag_name} -m "Add source sha"`
  puts("    Force annotation return code: #{$?}")
end
