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
RT_URL = "https://#{ENV['RT_URL']}"
RT_DIR = File.expand_path(ENV['RT_DIR'], CI_PATH)
ARCH = File.read(File.expand_path('build/arch.txt',CI_PATH))

############################################################
def main(ci_path, rt_url, rt_dir, arch, test_dir)
  debug = false
  g_ci = Git.open(ci_path)
  the_branch = determine_branch
  puts("  Current branch, #{the_branch}, obtained")
  the_ci_sha = g_ci.object("HEAD").sha
  puts("  SHA of main repository HEAD: #{the_ci_sha}")
  the_ci_msg = g_ci.object("HEAD").message
  puts("  Message of HEAD obtained:\n  #{the_ci_msg}")
  if get_tag
    results_branch = get_tag
  else
    results_branch = the_branch
  end
  puts("  - Results repo URL: #{rt_url}")
  puts("  - Results directory: #{rt_dir}")
  puts("  - Results branch: #{results_branch}")
  g_rt = robust_clone(rt_url, rt_dir, results_branch)
  if debug
    UTILS::git_status(g_rt.dir)
    UTILS::git_log(g_rt.dir)
  end
  puts("  RegressTest repo cloned to: #{g_rt.dir}")
  g_rt.config('user.name', "CI: #{arch}")
  g_rt.config('user.email', "ci@ci.org")
  `cd "#{g_rt.dir}"` # && git config --global credential.helper store`
  #File.write("#{ENV['HOME']}/.git-credentials", "https://$($env:PATOKEN):x-oauth-basic@github.com\n", mode: 'wb')
  mimic_source(g_rt, results_branch, the_ci_sha)
  puts("  Source mimicked")
  if debug
      UTILS::git_status(g_rt.dir)
      UTILS::git_log(g_rt.dir)
      UTILS::git_list_branches(g_rt.dir)
  end
end

puts("\n\n\nCloning test results repository.")

if not is_pull_request?
  main(CI_PATH, RT_URL, RT_DIR, ARCH, TEST_DIR)
else
  puts("Skipping clone and mimic due to being a pull request")
end

puts("Cloning test results repository completed!")
