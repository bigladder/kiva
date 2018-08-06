require_relative('lib')
require_relative('utils')

############################################################
# Required Inputs
# path to the Continuous Integration Source Code Git Repo
THIS_DIR = File.dirname(__FILE__)
CI_PATH = File.expand_path('..', THIS_DIR)
TEST_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
# regression testing repository URL
# PATOKEN = personal access token
RT_URL = "https://#{ENV['PATOKEN']}:x-oauth-basic@#{ENV['RT_URL']}"
RT_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
ARCH = File.read(File.expand_path('build/arch.txt',CI_PATH))

############################################################
def main(ci_path, rt_dir, arch, test_dir, rt_url)
  debug = false
  g_ci = Git.open(ci_path)
  the_branch = determine_branch
  puts("  Current branch, #{the_branch}, obtained")
  the_ci_sha = g_ci.object("HEAD").sha
  puts("  SHA of HEAD obtained: #{the_ci_sha}")
  the_ci_msg = g_ci.log[0].message # used to be "commit('HEAD')", but it creates bad message after GitHub merge
  puts("  Message of HEAD obtained:\n  #{the_ci_msg}")
  g_rt = Git.open(rt_dir)
  puts("  RegressTest repo opened at: #{g_rt.dir}")
  if debug
    UTILS::git_status(g_rt.dir)
  end
  FileUtils.cp(File.join('..', 'build', 'Testing', 'Temporary', 'LastTest.log'), File.join(rt_dir, arch))
  puts("  LastTest.log copied")
  puts("  Attempting to add any new files")
  chngs = UTILS.list_changes(g_rt.dir)
  puts("  Logging #{chngs['Deleted'].length} deleted files")
  chngs['Deleted'].each {|f| `cd #{g_rt.dir} && git rm --force #{f}`}
  g_rt.add(:all=>true)
  code = system("cd #{g_rt.dir} && git diff --quiet --cached --exit-code")
  if code
    puts("  No changes found")
  else
    puts("  Changes found. Committing...")
    the_commit = "#{the_ci_msg} [#{arch}]\n{:src-sha \"#{the_ci_sha}\" " +
      ":src-msg \"#{the_ci_msg}\" " +
      ":arch \"#{arch}\" " +
      ":src-branch \"#{the_branch}\"}"
    g_rt.commit(the_commit)
    puts("  Committed")
    puts("  Tagging...")
    tag_name = "src_#{the_ci_sha}_#{arch}"
    tag_exists = !(g_rt.tags.select {|t| t.name == tag_name}.empty?)
    if tag_exists
      puts("  Tag exists, retagging...")
      retag(g_rt.dir, tag_name, rt_url)
    else
      puts("  Adding tag #{tag_name}")
      g_rt.add_tag(tag_name, {a: true, m: "Add source sha"})
    end
    puts("  Tag added")
    puts("  Attempting to push/pull")
    robust_push_pull(g_rt, the_branch, the_commit, tag_name, rt_url,
                     our_dir=arch, to_be_deleted=chngs['Deleted'])
    puts("  Push/pull succeeded")
  end
  if debug
    UTILS::git_status(g_rt.dir)
    UTILS::git_log(g_rt.dir)
  end
end

puts("\n\n\nLogging results")

if not is_pull_request?
  main(CI_PATH, RT_DIR, ARCH, TEST_DIR, RT_URL)
else
  puts("Skipping logging results due to being a pull request")
end

puts("Results logged!")
