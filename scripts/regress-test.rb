require('git')
require('fileutils')
require('open3')

# String String String String -> Git::Base
# Robustly clones a source repo to the given target
# repo_url: String   = URL to repository to clone from 
# branch: String     = Commit SHA or Branch Name of a branch to switch to
# ref_commit: String = Commit SHA, Tag, or Branch Name to serve as the tag to
#                      branch off of
# local_path: String = File path to where the cloned repository should live on
#                      the local filesystem
def robust_clone(repo_url, local_path)
  g = nil
  if File.exist?(local_path)
    if File.exist?(File.join(local_path, '.git'))
      g = Git.open(local_path)
      r = g.remotes.select {|r| r.url == repo_url}[0]
      origin = g.remotes.select {|r| r.name == "origin"}[0]
      if r.nil? or r != origin
        g.remove_remote("origin") unless origin.nil?
        g.add_remote("origin", repo_url)
      end
    else
      if Dir[File.join(local_path, '*')].empty?
        g = Git.clone(repo_url, local_path)
      else
        raise "local_path \"#{local_path}\" exists but is NOT empty; please remove files and try again"
      end
    end
  else
    FileUtils.mkdir_p(local_path)
    g = Git.clone(repo_url, local_path)
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

def mimic_source(g, branch, src_commit)
  puts("Mimicing source repository")
  puts("- branch: #{branch}")
  if branch.nil? or branch.empty? or branch.include?("(") or branch.include?(" ")
    puts("detected nil/empty or malformed branch")
    puts("falling back to system reported branch name")
    if ! ENV['SRC_BRANCH'].nil?
      branch = ENV['SRC_BRANCH']
    elsif ! ENV['TRAVIS_BRANCH'].nil?
      puts("Tried to find environment variable SRC_BRANCH but it was nil...")
      puts("trying TRAVIS_BRANCH...")
      branch = ENV['TRAVIS_BRANCH']
      puts("TRAVIS_BRANCH found!")
    elsif ! ENV['APPVEYOR_REPO_BRANCH'].nil?
      puts("Tried to find environment variable SRC_BRANCH but it was nil...")
      puts("trying APPVEYOR_REPO_BRANCH...")
      branch = ENV['APPVEYOR_REPO_BRANCH']
      puts("APPVEYOR_REPO_BRANCH found!")
    end
  end
  puts("- branch: #{branch}")
  puts("- src_commit: #{src_commit}")
  ref_commit = find_commit_from_substr(g, src_commit)
  puts("- ref_commit = #{ref_commit}")
  if ref_commit.nil?
    puts("no relevant tag in repository")
  elsif not g.is_branch?(branch)
    puts("checking out upstream")
    g.checkout(ref_commit)
  end
  puts("creating and checking out new branch, #{branch}")
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
  w['stdout.log', stdout]
  w['stderr.log', stderr]
  w['exitcode.log', exitcode]
  c[:input_files].each do |f|
    FileUtils.rm(File.join(path, File.basename(f)))
  end
  FileUtils.rm(File.join(path, File.basename(a[:exe])))
  files_to_add
end

# robust pull/push
def robust_push_pull(g, branch)
  puts("Starting robust_pull_push!")
  begin
    puts("Attempting to pull")
    g.pull('origin', branch) if g.is_remote_branch?(branch)
    puts("Pull attempt succeeded")
  rescue => e
    puts("Trying to fix suspected auto-merge conflict")
    puts("Error: #{e.message}")
    puts("Git repo directory: #{g.dir}")
    # OK, we probably have a merge conflict
    # we need to:
    # 1. find which files are in conflict
    # 2. `git checkout --ours #{filename}` for each file in conflict
    # 3. re-add the files
    # 4. re-commit the files with some sort of commit message

    # `git diff --name-only --diff-filter=U`
    puts("Asking git for list of conflicted files")
    files_in_conflict = `cd #{g.dir} && git diff --name-only --diff-filter=U`
    puts("files in conflict: #{files_in_conflict.split.inspect}")
    files_in_conflict.split.each do |fname|
      puts("checking out and re-adding #{fname}")
      `cd #{g.dir} && git checkout --ours #{fname}`
      g.add(fname)
    end
    puts("(Re-)Committing...")
    g.commit("Commit to fix autoconflict")
    puts("Done!")
  end
  puts("Pushing to origin!")
  g.push('origin', branch, {:tags=>true})
  puts("Done robust_pull_push!")
end

# From CI
# - what branch is CI on
# - commit-sha
# - commit-message
# - 
# NOTE: used tags in regresstest repo to identify which commit in source a regress-test corresponds to

############################################################
# Required Inputs
# path to the Continuous Integration Source Code Git Repo
THIS_DIR = File.dirname(__FILE__)
CI_PATH = File.expand_path('..', THIS_DIR)
TEST_DIR = File.expand_path(ENV['TEST_DIR'], CI_PATH)
# regression testing repository URL
# PATOKEN = personal access token
RT_URL = "https://#{ENV['PATOKEN']}@#{ENV['RT_URL']}"
RT_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
ARCH = ENV['BUILD_ARCHITECTURE']
#APP = {
#  exe: "kiva.rb",
#  cmd: "ruby kiva.rb"
#}
#CASES = [
#  {
#    id: "case-1",
#    app: APP,
#    input_files: ["tests/case-1/input.yaml", "tests/chicago.epw"],
#    output_files: ["out.csv"],
#    args: ["input.yaml", "chicago.epw", "out.csv"]
#  },
#  {
#    id: "case-2",
#    app: APP,
#    input_files: ["tests/case-2/input.yaml", "tests/chicago.epw"],
#    output_files: ["out.csv"],
#    args: ["input.yaml", "chicago.epw", "out.csv"],
#  }
#]

############################################################
def main(ci_path, rt_url, rt_dir, arch, test_dir)
  puts("Starting main")
  puts("Opening git on CI_PATH")
  g_ci = Git.open(ci_path)
  puts("Opened...")
  puts("Getting current branch...")
  the_branch = g_ci.current_branch # current branch
  if the_branch.include?("(") or the_branch.include?(" ")
    the_branch = ENV['TRAVIS_BRANCH']
  end
  puts("Current branch, #{the_branch}, obtained")
  puts("Getting SHA of HEAD...")
  the_ci_sha = g_ci.object("HEAD").sha
  puts("Sha of HEAD obtained")
  puts("Getting Message of HEAD")
  the_ci_msg = g_ci.object("HEAD").message
  puts("Message of HEAD obtained")
  puts("Attempting to CLONE the RegressTest repo")
  g_rt = robust_clone(rt_url, rt_dir)
  puts("RegressTest repo cloned")
  puts("- RegressTest repo directory: #{g_rt.dir}")
  puts("Setting Git username and email")
  g_rt.config('user.name', "CI: #{arch}")
  g_rt.config('user.email', "ci@ci.org")
  puts("Git username and email set")
  puts("Attempting to mimic source")
  mimic_source(g_rt, the_branch, the_ci_sha)
  puts("Source mimicked")
  # Run Kiva
  puts("Copying case files to repo")
  puts("test_dir = #{test_dir}")
  puts("test_dir exist? #{File.exist?(test_dir)}")
  #f = lambda do |p|
  #  puts("contents of #{p}: #{Dir[File.join(p, '*')].sort}")
  #end
  #f["/home/travis/build/michael-okeefe/kiva/build/test"]
  #f["/home/travis/build/michael-okeefe/kiva/build/Testing"]
  FileUtils.cp_r(File.join(test_dir, '.'), rt_dir)
  puts("Case files copied")
  #files_to_add = []
  #cases.each do |c|
  #  a = c[:app]
  #  puts("running case #{c[:id]}")
  #  files_to_add += run_case(ci_path, rt_dir, arch, a, c)
  #end
  # Copy files to copy to the checked out repository 
  # Set files_to_add
  #puts("Adding #{files_to_add.length} files")
  puts("Attempting to see if there are any changes")
  code = `cd #{g_rt.dir} && git diff --exit-code`
  if code.empty?
    puts("No changes found")
    puts("Output of git status is:")
    puts(`cd #{g_rt.dir} && git status`)
  else
    puts("Changes found")
    puts("Attempting to add all files")
    g_rt.add(:all=>true)
    puts("All files added")
    puts("Committing...")
    g_rt.commit(
      "#{the_ci_msg} [#{arch}]\n{:src-sha \"#{the_ci_sha}\" " +
      ":src-msg \"#{the_ci_msg}\" " +
      ":arch \"#{arch}\" " + 
      ":src-branch \"#{the_branch}\"}"
    )
    puts("Committed")
    puts("Tagging...")
    tag_name = "src_#{the_ci_sha}"
    tag_exists = ! g_rt.tags.select {|t| t.name == tag_name}.empty?
    if tag_exists
      puts("Tag, #{tag_name}, exists")
      # delete tag on remote
      `cd #{g_rt.dir} && git push origin :refs/tags/#{tag_name}`
      # force annotate the tag again
      `cd #{g_rt.dir} && git tag -fa #{tag_name} -m "Add source sha"`
      # push to origin will occur in a bit
    else
      puts("Adding tag #{tag_name}")
      g_rt.add_tag(tag_name, {a: true, m: "Add source sha"})
    end
    puts("Tag added")
    puts("Attempting to push/pull")
    robust_push_pull(g_rt, the_branch)
    puts("push/pull succeeded")
  end
end
main(CI_PATH, RT_URL, RT_DIR, ARCH, TEST_DIR)
puts("Done!")
