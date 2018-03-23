
# only prints if str is nonempty
function printiferror(log,str)
	str = strip(str)
	isempty(str) || println(log, "ERROR: ", str)
end

# only prints if str is nonempty
function printifwarning(log,str)
	str = strip(str)
	isempty(str) || println(log, "WARNING: ", str)
end

# only prints if str is nonempty
function printifinfo(log,str)
	str = strip(str)
	isempty(str) || println(log, str)
end




# annoying helper functions that are needed since I cannot redirect stdout/err to IOBuffers
# and working with Pipes directly is difficult as they might block, waiting for input, 
# and it's hard to know when that will happen.
# Should be replaced by a better solution using IOBuffers or Pipes directly when Julia has better support.

# TODO: fix using redirect_stdout/redirect_stderr ???


type DiskBuffer
	filename::String
	stream#::IOStream
end
DiskBuffer() = DiskBuffer(tempname(),Void) 

openbuf(db::DiskBuffer) = (db.stream = open(db.filename, "w"); db.stream::IOStream)
function closebuf(db::DiskBuffer)
	close(db.stream)
	str = readstring(db.filename)
	rm(db.filename) # delete file
	str
end


