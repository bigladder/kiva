require_relative('lib')

############################################################
# Required Inputs
# path to the Continuous Integration Source Code Git Repo
THIS_DIR = File.dirname(__FILE__)
CI_PATH = File.expand_path('..', THIS_DIR)
TEST_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
# regression testing repository URL
# PATOKEN = personal access token
RT_URL = "https://#{ENV['PATOKEN']}@#{ENV['RT_URL']}"
RT_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
ARCH = ENV['BUILD_ARCHITECTURE']

############################################################
def main(ci_path, rt_url, rt_dir, arch, test_dir)
  puts("Starting main")
  puts("Opening git on CI_PATH")
  g_ci = Git.open(ci_path)
  puts("Opened...")
  puts("Getting current branch...")
  the_branch = determine_branch
  puts("Current branch, #{the_branch}, obtained")
  puts("Getting SHA of HEAD...")
  the_ci_sha = g_ci.object("HEAD").sha
  puts("Sha of HEAD obtained: #{the_ci_sha}")
  puts("Getting Message of HEAD")
  the_ci_msg = g_ci.object("HEAD").message
  puts("Message of HEAD obtained:\n#{the_ci_msg}")
  puts("Attempting to OPEN the RegressTest repo")
  g_rt = Git.open(rt_dir)
  puts("RegressTest repo opened")
  puts("- RegressTest repo directory: #{g_rt.dir}")
  puts("Attempting to add any new files")
  g_rt.add(:all=>true)
  puts("All files added")
  puts("Attempting to see if there are any changes cached after 'add --all'")
  code = system("cd #{g_rt.dir} && git diff --quiet --cached --exit-code")
  if code
    puts("No changes found")
    puts("- code is: #{code}")
  else
    puts("Changes found")
    puts("Committing...")
    the_commit = "#{the_ci_msg} [#{arch}]\n{:src-sha \"#{the_ci_sha}\" " +
      ":src-msg \"#{the_ci_msg}\" " +
      ":arch \"#{arch}\" " +
      ":src-branch \"#{the_branch}\"}"
    g_rt.commit(the_commit)
    puts("Committed")
    puts("Tagging...")
    tag_name = "src_#{the_ci_sha}_#{arch}"
    tag_exists = ! g_rt.tags.select {|t| t.name == tag_name}.empty?
    if tag_exists
      puts("Tag, #{tag_name}, exists")
      # delete tag on remote
      `cd #{g_rt.dir} && git push --quiet origin :refs/tags/#{tag_name}`
      # force annotate the tag again
      `cd #{g_rt.dir} && git tag -fa #{tag_name} -m "Add source sha"`
      # push to origin will occur in a bit
    else
      puts("Adding tag #{tag_name}")
      g_rt.add_tag(tag_name, {a: true, m: "Add source sha"})
    end
    puts("Tag added")
    puts("Attempting to push/pull")
    robust_push_pull(g_rt, the_branch, the_commit)
    puts("push/pull succeeded")
  end
end

puts("scripts/log-results.rb Start!")

main(CI_PATH, RT_URL, RT_DIR, ARCH, TEST_DIR)

puts("scripts/log-results.rb Done!")
