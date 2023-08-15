# Applies the boost bcp utility to the src folder of a repo.
# This will replace the current boost folder (already in vendor folder of the repo)
# with a new boost folder that should contain only the need portions of boost,
# and which may be retrieved from a different version of boost.
# The previous boost folder will be stashed in the same location, with
# a different name, if necessary, and can be deleted after the transition is confirmed.
#
# Notes:
# bcp is not able to trace dependencies through macros, resulting in some 'CAUTION' statements.
# In the case of kiva, these do not appear to be consequential.
# -----------------------------------
# bcp can be installed/updated on Mac using homebrew, prior to running this script:
#   brew install boost-bcp
#   brew link boost-bcp
#
# Example of use on Mac:
# 	python3 update-boost.py bcp /Users/USERNAME/Documents/Development/boost_1_81_0 \
#		/Users/USERNAME/Documents/GitHub/kiva boost-1.77.0 boost-1.81.0
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
#	(A large number of messages is displayed.)
#
# The bcp Unix executable is in:
#	boost_X_XX_X/dist/bin/bcp
# and also in (for some reason):
#	bin.v2/..
#
# Example of use on Ubuntu (or Mac):
# 	python3 update-boost.py /Users/USERNAME/Documents/Development/boost_1_81_0/dist/bin/bcp \
# 		/Users/USERNAME/Documents/Development/boost_1_81_0 \
# 		/Users/USERNAME/Documents/GitHub/kiva boost-1.77.0 boost-1.81.0
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
#	(Large number of messages displayed.)
#
# bcp.exe is in:
#	boost_X_XX_X/dist/bin/bcp
# and also in (for some reason):
#	bin.v2/..
#
# Example of use on Windows (e.g., PowerShell):
# 	python3 update-boost.py \
#		"C:\Users\USERNAME\Documents\Development\boost_1_81_0\dist\bin\bcp.exe" \
#		"C:\Users\USERNAME\Documents\Development\boost_1_81_0" \
# 		"C:\Users\USERNAME\Documents\kiva" boost-1.77.0 boost-1.81.0

import os
import sys
import subprocess
import shutil

# replace a word in a file
def find_and_replace(filename, word, replacement):
	text = ''
	with open(filename, 'r') as f:
		text = f.read()
	text = text.replace(word, replacement)
	with open(filename, 'w') as f:
		f.write(text)			

# Find a unique folder name and rename the folder
def check_rename_folder(orig_path):
	i = 0
	while os.path.exists(orig_path + "_" + str(i)):
		i = i + 1
	mod_path = orig_path + "_" + str(i)

	# rename the folder
	print("Renaming folder: " + orig_path + " to " + mod_path)
	os.rename(orig_path, mod_path)
	return mod_path

# source_boost_path: path to folder containing new, full boost version
# repo_path: path to repo
# prev_repo_boost_folder_name: name of current boost folder (within repo)
# dest_repo_boost_folder_name: new name for updated, minimized boost folder (within repo)
def updateBoost(bcp_command, source_boost_path, repo_path, prev_repo_boost_folder_name, dest_repo_boost_folder_name):

	# path to folder and list of sub-folders containing source code that uses boost
	repo_src_path = os.path.join(repo_path, "src")
	
	# generate full paths
	prev_repo_boost_path = os.path.join(repo_path, "vendor", prev_repo_boost_folder_name)
	dest_repo_boost_path = os.path.join(repo_path, "vendor", dest_repo_boost_folder_name)

	# name with which to stash current boost folder (which can then be deleted)	
	stash_boost_path = prev_repo_boost_path

	# modify stash folder name if same as dest
	if prev_repo_boost_path == dest_repo_boost_path:
		stash_boost_path = check_rename_folder(prev_repo_boost_path)

	# if dest folder name is used, rename folder that currently has that name
	if os.path.exists(dest_repo_boost_path):
		check_rename_folder(dest_repo_boost_path)

	# make new dest repo boost folder
	os.mkdir(dest_repo_boost_path )
	
	# run bcp utility
	bcp_list = [bcp_command, "--scan"]
	bcp_list.append("--boost=" + source_boost_path)
	for repo_src_path, dirs, files in os.walk(repo_src_path, topdown = True):
		for filename in files:
			full_filename = os.path.join(repo_src_path, filename)
			bcp_list.append(full_filename)
	
	bcp_list.append(dest_repo_boost_path)
	
	result = subprocess.run(bcp_list, stdout = subprocess.PIPE, text = True)
	print(result.stdout)
	
	# copy CMakeLists.txt from stashed boost folder
	full_filename = os.path.join(stash_boost_path, "CMakeLists.txt")
	if os.path.exists(full_filename):
		print("Copying file: " + full_filename)
		shutil.copyfile(full_filename, os.path.join(dest_repo_boost_path, "CMakeLists.txt"))
	else:
		print(full_filename + " does not exist")
	
	# copy CMakeLists.txt from stashed boost program_options folder
	full_filename = os.path.join(stash_boost_path, "libs", "program_options", "CMakeLists.txt")
	if os.path.exists(full_filename):
		dest_folder = os.path.join(dest_repo_boost_path, "libs", "program_options")
		if os.path.exists(dest_folder):
			print("Copying file: " + full_filename)
			shutil.copyfile(full_filename, os.path.join(dest_folder, "CMakeLists.txt"))
		else:
			print(dest_folder + " does not exist")
	else:
		print(full_filename + " does not exist")
		
	# change the boost folder name in the CMakeLists.txt file in the vendor folder
	full_filename = os.path.join(repo_path, "vendor", "CMakeLists.txt")
	if os.path.exists(full_filename):
		print("Changing boost reference in " + full_filename)
		find_and_replace(full_filename, prev_repo_boost_folder_name, dest_repo_boost_folder_name)
	else:
		print(full_filename + " does not exist")	

	# copy license
	file_list = os.listdir(source_boost_path)
	for filename in file_list:
		if "LICENSE" in filename:
			license_name = filename
			full_filename = os.path.join(source_boost_path, license_name)
			print("Copying file: " + full_filename)
			shutil.copyfile(full_filename, os.path.join(dest_repo_boost_path, license_name))	
			
#  main
n_args = len(sys.argv) - 1

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
