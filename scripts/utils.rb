
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
end
