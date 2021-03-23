require 'fileutils'

FileUtils.rm_rf "./xml"
FileUtils.mkdir "./xml"
puts "Start"
Dir.glob("./xml-*/*").each do |from|
  puts "Copy: #{from}"
  to = "./xml/#{File.basename(from)}"
  if File.exists?(to) then
    system("xml-cat #{from} #{to} > tmp");
    FileUtils.mv "tmp", to 
  else
    FileUtils.cp from, to 
  end
end
puts "Done"
