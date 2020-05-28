using Weave
pwd = "/Users/christianhilscher/Desktop/Assignments/TSA-Assignments/Assignment10"
out = "/Users/christianhilscher/Desktop/Assignments/TSA-Assignments/Assignment10/output"
# add depencies for the example

filename = normpath(pwd, "Assignment10.jmd")
#weave(filename1, out_path = :pwd,doctype = "md2pdf")
weave(filename, out_path = out,doctype = "md2pdf")
