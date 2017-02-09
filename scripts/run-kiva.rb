require('open3')
require('fileutils')
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
  FileUtils.rm_rf(output_dir) if File.exist?(output_dir)
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
    "\"#{run_pat}#{File.basename(exe_path)}\" \"#{in_file}\" \"#{weather_file}\" \"#{output_path}\""
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
  return exitcode.success?
end

KIVA_PATH = File.expand_path(ARGV[0])
KIVA_DIR = File.dirname(KIVA_PATH)
INPUT_FILE = File.expand_path(ARGV[1])
WEATHER_FILE = File.expand_path(ARGV[2])
OUTPUT_FILE  = File.expand_path(ARGV[3])
OUTPUT_DIR = File.dirname(OUTPUT_FILE)
puts("Running with:")
puts("  kiva path  = #{KIVA_PATH}")
puts("  input file = #{INPUT_FILE}")
puts("  weather file = #{WEATHER_FILE}")
puts("  output file  = #{OUTPUT_FILE}")
success = run_case(KIVA_PATH, INPUT_FILE, WEATHER_FILE, OUTPUT_FILE)
f = lambda do |dir|
  puts("Evaluating contents of #{dir}")
  if File.exists?(dir)
    puts("- contents:\n  #{Dir[File.join(dir, '*')]}")
  else
    puts("- #{dir} doesn't exist...")
  end
end
f[OUTPUT_DIR]
f[KIVA_DIR]
puts("run-kiva.rb completed!")

exit(success)
