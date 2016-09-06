require('open3')
require('time')

puts("Running Kiva!!!")
THIS_DIR = File.expand_path(File.dirname(__FILE__))

def operating_system
  if (/darwin/ =~ RUBY_PLATFORM)
    :mac
  elsif (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM)
    :windows
  else
    :unix
  end
end

def run_case(exe_path, in_file, weather_file, output_path)
  output_dir = File.dirname(output_path)
  puts("  ... output directory=#{output_dir}")
  puts("  ... exe path        =#{exe_path}")
  FileUtils.mkdir_p(output_dir) unless File.exist?(output_dir)
  FileUtils.cp(in_file, output_dir)
  FileUtils.cp(weather_file, output_dir)
  run_pat = if operating_system == :windows
              ""
            else
              "./"
            end
  cmd = [
    "cd #{File.dirname(exe_path)}",
    "#{run_pat}#{File.basename(exe_path)} #{in_file} #{weather_file} #{output_path}"
  ].join(" && ")
  puts("  ... cmd = #{cmd}")
  t_start = Time.now
  stdout, stderr, exitcode = Open3.capture3(cmd)
  t_end = Time.now
  w = lambda do |name, data|
    p = File.join(output_dir, name)
    File.write(p, data)
  end
  w['stdout.log', stdout]
  w['stderr.log', stderr]
  w['exitcode.log', exitcode]
  w['timings.log', (t_end - t_start).to_s + " seconds"]
  puts("run_case complete!")
end

KIVA_PATH = File.expand_path(ARGV[0])
INPUT_FILE = File.expand_path(ARGV[1])
WEATHER_FILE = File.expand_path(ARGV[2])
OUTPUT_FILE  = File.expand_path(ARGV[3])
puts("Running with:")
puts("  kiva path  = #{KIVA_PATH}")
puts("  input file = #{INPUT_FILE}")
puts("  weather file = #{WEATHER_FILE}")
puts("  output file  = #{OUTPUT_FILE}")
run_case(KIVA_PATH, INPUT_FILE, WEATHER_FILE, OUTPUT_FILE)
