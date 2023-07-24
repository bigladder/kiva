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
# -----------------------------------
# bcp can be installed/updated on Mac using homebrew, prior to running this script:
#   brew install boost-bcp
#   brew link boost-bcp
#
# Example of use on Mac:
# python3 update-boost.py bcp /Users/USERNAME/Documents/Development/boost_1_81_0 \
#	/Users/USERNAME/Documents/GitHub/kiva boost-1.77.0 boost-1.81.0
# -----------------------------------
# bcp can be installed on Ubuntu (or Mac) as follows:
# 1) Download boost version (boost_X_XX_X) at www.boost.org/users/download/
#
# 2) In boost_X_XX_X folder, on command line, enter:
# 	./bootstrap.sh
#	(Lots of warning, etc., displayed.)
#	(b2 Unix executable appears in boost_X_XX_X folder.)
#
# 3) Enter:
#	./b2 tools/bcp
#	(A large number of messages is displayed)
#
# The bcp Unix executable is in:
#	boost_X_XX_X/dist/bin/bcp
# and also in (for some reason):
#	bin.v2/..
#
# Example of use on Ubuntu (or Mac):
# python3 update-boost.py /Users/USERNAME/Documents/Development/boost_1_81_0/dist/bin/bcp \
# /Users/USERNAME/Documents/Development/boost_1_81_0 \
# /Users/USERNAME/Documents/GitHub/kiva boost-1.77.0 boost-1.81.0
#
# -------------------------------------
# bcp can be installed on Windows as follows:
# 1) Download boost version (boost_X_XX_X) at www.boost.org/users/download/
#
# 2) In a DOS Shell, go to boost_X_XX_X folder, and enter:
#	bootstrap.bat
#	(Lots of warning, etc., displayed.)
#	(b2.exe appears in boost_X_XX_X folder.)
#
# 3) Enter:
#	b2 tools/bcp
#	(large number of messages displayed)
#
# bcp.exe is in:
#	boost_X_XX_X/dist/bin/bcp
# and also in (for some reason):
#	bin.v2/..
#
# Example of use on Windows (e.g., PowerShell):
# python3 update-boost.py \
#	"C:\Users\USERNAME\Documents\Development\boost_1_81_0\dist\bin\bcp.exe" \
#	"C:\Users\USERNAME\Documents\Development\boost_1_81_0" \
# 	"C:\Users\USERNAME\Documents\kiva" boost-1.82.0 boost-1.81.0

import os
import sys
import re
import platform
import subprocess
import shlex

# replace a word in a file
def find_and_replace(file, word, replacement):
    text = ''
    with open(file, 'r') as f:
        text = f.read()
    text = text.replace(word, replacement)
    with open(file, 'w') as f:
      f.write(text)

def add_files(file_list, folder_path):
	path_items = os.listdir(foler_path)
			

# source_boost_path: path to folder containing new, full boost version
# repo_path: path to repo
# prev_repo_boost_folder_name: name of current boost folder (within repo)
# dest_repo_boost_folder_name: new name for updated, minimized boost folder (within repo)
def updateBoost(bcp_command, source_boost_path, repo_path, prev_repo_boost_folder_name, dest_repo_boost_folder_name):

	# path to folder and list of subfolders containing source code that uses boost
	repo_src_path = repo_path + "/src"
	
	# name with which to stash current boost folder (which can then be deleted)	
	stash_repo_boost_folder_name = prev_repo_boost_folder_name
	if stash_repo_boost_folder_name == dest_repo_boost_folder_name:
		stash_repo_boost_folder_name += "_old"

	# rename old boost folder
	prev_repo_boost_path = repo_path + "/vendor/" + prev_repo_boost_folder_name
	dest_repo_boost_path = repo_path + "/vendor/" + dest_repo_boost_folder_name
	stash_boost_path = repo_path + "/vendor/" + stash_repo_boost_folder_name

	# append a number if stash folder exists
	if stash_repo_boost_folder_name != prev_repo_boost_folder_name:
		if os.path.exists(stash_boost_path):
			i = 0
			while os.path.exists(stash_boost_path + "_" + str(i)):
				i = i + 1
			stash_boost_path = stash_boost_path + "_" + str(i)

		# rename the current repo boost folder
		os.system("mv \"" + prev_repo_boost_path + "\" \"" + stash_boost_path + "\"")

	# make new repo boost folder
	if not(os.path.exists(dest_repo_boost_path)):
		os.system("mkdir \"" + dest_repo_boost_path + "\"")

	# run bcp utility
	file_list = ""
	for repo_src_path, dirs, files in os.walk(repo_src_path, topdown = True):
		for filename in files:
			file_list = file_list + " \"" + os.path.join(repo_src_path, filename) + "\""
			
	cwd = os.getcwd()
	if platform.system() == "Windows":
		bcp_dir = os.path.dirname(bcp_command)
		bcp_command = os.path.basename(bcp_command)
		if os.path.exists(bcp_dir):
			os.chdir(bcp_dir)
		else:
			quit()

	full_bcp_command = bcp_command + " --scan"
	full_bcp_command = full_bcp_command + " --boost="
	full_bcp_command = full_bcp_command + "\"" + source_boost_path + "\""

	full_bcp_command = full_bcp_command + file_list
	full_bcp_command = full_bcp_command + " \"" + dest_repo_boost_path +"\""

	#print(full_bcp_command)
	os.system(full_bcp_command)

	if platform.system() == "Windows":
		os.chdir(cwd)
	
	copy_command = "cp"
	if platform.system() == "Windows":
		copy_command = "copy"
		stash_boost_path  = stash_boost_path.replace("/", "\\")
 
	# copy CMakeLists.txt from stashed boost folder and subfolder(s)
	if os.path.exists(stash_boost_path + "/CMakeLists.txt"):
		print("Copying  \"" + stash_boost_path + "/CMakeLists.txt\"")
		os.system(copy_command + " \"" + stash_boost_path + "/CMakeLists.txt\" \"" + dest_repo_boost_path + "\"")
	else:
		print("\"" + stash_boost_path + "/CMakeLists.txt\""+ " does not exist")

	if os.path.exists(stash_boost_path + "/libs/program_options/CMakeLists.txt"):
		if os.path.exists(dest_repo_boost_path + "/libs/program_options"):
			print("Copying  \"" + source_boost_path + "/libs/program_options/CMakeLists.txt\"")
			os.system(copy_command + " \"" + stash_boost_path + "/libs/program_options/CMakeLists.txt\" \"" + dest_repo_boost_path + "/libs/program_options\"")
		else:
			print("\"" + source_boost_path + "/libs/program_options/CMakeLists.txt\"" + " does not exist")

	# change the boost folder name in the CMakeLists.txt file in the vendor folder
	if os.path.exists(repo_path + "/vendor/CMakeLists.txt"):
		print("Changing boost reference in  \"" + repo_path + "/vendor/CMakeLists.txt\"")
		find_and_replace(repo_path + "/vendor/CMakeLists.txt", prev_repo_boost_folder_name, dest_repo_boost_folder_name)

	# copy license
	if os.path.exists(source_boost_path + "/LICENSE*.*"):
		print("Copying  \"" + source_boost_path + "/LICENSE*.*\"")
		if platform.system() == "Windows":
			os.system(copy_command + " \"" + source_boost_path + "/LICENSE*.*\" \"" + dest_repo_boost_path + "\"")
		else:
			os.system(copy_command + " " + source_boost_path + "/LICENSE*.* " + dest_repo_boost_path)
	
#  main
n_args = len(sys.argv) - 1

for i in range(1, n_args + 1):
	print(sys.argv[i])
	
if n_args == 5:
	bcp_command = sys.argv[1]
	source_boost_path = sys.argv[2]
	repo_path = sys.argv[3]
	prev_repo_boost_folder_name = sys.argv[4]
	dest_repo_boost_folder_name = sys.argv[5]

	updateBoost(bcp_command, source_boost_path, repo_path, prev_repo_boost_folder_name, dest_repo_boost_folder_name)

else:
	print('arguments:')
	print('1. command to invoke bcp utility')
	print('2. path to root of boost library to install')
	print('3. path to root of repo')
	print('4. name of previous boost folder to replace (will be renamed, if needed)')
	print('5. name of new boost folder to install')
