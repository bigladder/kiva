
module UTILS
  def self.universal_newline(s)
    s.encode(s.encoding, universal_newline: true)
  end

  def self.list_changes(dir, printout=true)
    out = {}
    [['A', 'Added'], ['C', 'Copied'], ['D', 'Deleted'], ['M', 'Modified'],
     ['R', 'Renamed'], ['T', 'Type Changed'], ['U', 'Unmerged'],
     ['X', 'Unknown'], ['B', 'Broken']].each do |key, tag|
       files = `cd #{dir} && git diff --name-only --diff-filter=#{key}`
       out[tag] = files.split
       if printout
         puts(tag + ':')
         puts(files)
         puts('')
       end
     end
     out
  end

  def self.git_status(dir)
    puts("Calling git status")
    puts(`cd #{dir} && git status`)
  end

  def self.git_log(dir)
    puts("Calling git log")
    puts(`cd #{dir} && git log --max-count=10 --pretty=format:"%h - {%p} [%an]: %s"`)
  end

  def self.git_list_branches(dir)
    puts("Detailed List of Local and Remote Branches:")
    puts(`cd #{dir} && git branch -a -v`)
  end

  def self.git_remote_branch_exists?(remote_url, branch)
    `git ls-remote --exit-code --heads #{remote_url} #{branch}`
    $?.exitstatus == 0
  end
end
