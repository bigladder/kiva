# Applies the boost bcp utility on specified portion of a repo
# This will replace the current boost folder (already in the repo)
# with a new boost folder that should contain only the need portions of boost,
# and which may be retrieved from a different version of boost.
# The previous version of boost will be stashed in the same location, with
# a different name, and can be deleted after the transition is confirmed.
#
# Notes:
# 1. The contents of the reduced boost folder will only be expanded (not reduced)
# if bcp is applied to different source-code selections.
#
# 2. bcp is not able to trace dependencies through macros, resulting in some 'CAUTION' statements.
# In the case of kiva, these do not appear to be consequential.
#
# bcp can be installed/updated using homebrew, prior to running this script:
#   brew install boost-bcp
#   brew link boost-bcp

import os

docs_path = os.path.expanduser("~/Documents")

# specifies command to execute boost utility
bcp_command = "bcp"

# path to folder containing new, full boost version
source_boost_path = docs_path + "/Development/boost_1_82_0"

# path to repo
repo_path = docs_path + "/GitHub/kiva"

# path to folder and list of subfolders containing source code that uses boost
repo_src_path=repo_path+"/src"
repo_src_folder_names=["kiva", "libgroundplot", "libkiva"]

# name of current boost folder (within repo)
prev_repo_boost_folder_name = "boost-1.77.0"

# new name for updated, minimized boost folder (within repo)
dest_repo_boost_folder_name = "boost-1.82.0"

# new name with which to stash current boost folder (which can then be deleted)
stash_repo_boost_folder_name = "old_boost"

# rename old boost folder
prev_repo_boost_path = repo_path + "/vendor/" + prev_repo_boost_folder_name
dest_repo_boost_path = repo_path + "/vendor/" + dest_repo_boost_folder_name
stash_boost_path = repo_path + "/vendor/" + stash_repo_boost_folder_name

# replace a word in a file
def find_and_replace(file, word, replacement):
    text = ''
    with open(file, 'r') as f:
        text = f.read()
    text = text.replace(word, replacement)
    with open(file, 'w') as f:
      f.write(text)

# append a number if stash folder exists
if os.path.exists(stash_boost_path):
    i = 0;
    while os.path.exists(stash_boost_path + "_" + str(i)):
        i = i + 1
    stash_boost_path = stash_boost_path + "_" + str(i)

# rename the current repo boost folder
os.system("mv \"" + prev_repo_boost_path + "\" \"" + stash_boost_path + "\"")

# make new repo boost folder
os.system("mkdir \"" + dest_repo_boost_path + "\"")

# run bcp utility
for src_folder in repo_src_folder_names:
    src_path = repo_src_path + "/" + src_folder + "/*.*"
    os.system(bcp_command + " --scan --boost=" + source_boost_path + " " + src_path+ " " + dest_repo_boost_path)

# copy CMakeLists.txt from stashed boost folder and subfolder(s)
os.system("cp " + stash_boost_path + "/CMakeLists.txt " + dest_repo_boost_path)
os.system("cp " + stash_boost_path + "/libs/program_options/CMakeLists.txt " + dest_repo_boost_path + "/libs/program_options")

# change the boost folder name in the CMakeLists.txt file in the vendor folder
find_and_replace(repo_path + "/vendor/CMakeLists.txt",prev_repo_boost_folder_name,dest_repo_boost_folder_name)

# copy license
os.system("cp " + source_boost_path + "/LICENSE*.* " + dest_repo_boost_path)
