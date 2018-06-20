require_relative('lib')

############################################################
# Required Inputs
# path to the Continuous Integration Source Code Git Repo
THIS_DIR = File.dirname(__FILE__)
CI_PATH = File.expand_path('..', THIS_DIR)
TEST_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
# regression testing repository URL
# PATOKEN = personal access token
RT_URL = "https://#{ENV['RT_URL']}"
RT_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
ARCH = File.read(File.expand_path('build/arch.txt',CI_PATH))

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
  puts("Attempting to CLONE the RegressTest repo")
  g_rt = robust_clone(rt_url, rt_dir)
  puts("RegressTest repo cloned")
  puts("- RegressTest repo directory: #{g_rt.dir}")
  puts("Setting Git username and email")
  g_rt.config('user.name', "CI: #{arch}")
  g_rt.config('user.email', "ci@ci.org")
  puts("Changing directory to mimic source: #{g_rt.dir}")
  `cd "#{g_rt.dir}"` # && git config --global credential.helper store`
  #File.write("#{ENV['HOME']}/.git-credentials", "https://$($env:PATOKEN):x-oauth-basic@github.com\n", mode: 'wb')
  puts("Git username and email set")
  puts("Attempting to mimic source")
  mimic_source(g_rt, the_branch, the_ci_sha)
  puts("Source mimicked")
end

puts("scripts/clone-and-mimic.rb Start!")

if not is_pull_request?
  main(CI_PATH, RT_URL, RT_DIR, ARCH, TEST_DIR)
else
  puts("Skipping clone and mimic due to being a pull request")
end

puts("scripts/clone-and-mimic.rb Done!")
