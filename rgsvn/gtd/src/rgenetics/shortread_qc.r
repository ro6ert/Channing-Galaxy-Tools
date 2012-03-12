library('ShortRead')
ca = commandArgs(T) 
qa_data = qa('/tmp', ca[1], ca[2])
report(qa_data, dest=ca[3])
