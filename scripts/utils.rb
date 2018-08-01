
module UTILS
  def self.universal_newline(s)
    s.encode(s.encoding, universal_newline: true)
  end
end
