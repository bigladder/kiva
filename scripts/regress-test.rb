require('git')
require('fileutils')
require('open3')
require('pry')

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
  branch = ENV['TRAVIS_BRANCH'] if branch.include?("(") or branch.include?(" ")
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
  g.pull('origin', branch) if g.is_remote_branch?(branch)
  g.push('origin', branch, {:tags=>true})
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
TEST_DIR = File.expand_path('../build/test/results')
# regression testing repository URL
# PATOKEN = personal access token
RT_URL = "https://#{ENV['PATOKEN']}@github.com/michael-okeefe/test.git"
RT_DIR = File.expand_path("regress", THIS_DIR)
ARCH = "#{ENV['TRAVIS_OS_NAME']}-#{ENV['PROCESSOR_ARCHITECTURE']}-#{ENV['COMPILER']}"
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
  puts("Setting Git username and email")
  g_rt.config('user.name', "TRAVIS CI: #{arch}")
  g_rt.config('user.email', "travis@travis.org")
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
  tag_name = "src_#{the_ci_sha}__arch_#{arch}"
  tag_exists = ! g_rt.tags.select {|t| t.name == tag_name}.empty?
  if tag_exists
    puts("Tag, #{tag_name}, exists")
  else
    puts("Adding tag #{tag_name}")
  end
  g_rt.add_tag(tag_name, {a: true, m: "Add source sha"}) unless tag_exists
  puts("Tag added")
  puts("Attempting to push/pull")
  robust_push_pull(g_rt, the_branch)
  puts("push/pull succeeded")
end
main(CI_PATH, RT_URL, RT_DIR, ARCH, TEST_DIR)
puts("Done!")
