f = "/home/yosukefujii/Desktop/par.txt"
outwd = "/home/yosukefujii/Desktop/par/"
import os
g = open(f, "rU")

pars = ["gcv", "gcf", "rhov", "rhof", "lambda", "mcn", "area", "vol", "rho"]
out = 0
w0 = map(lambda x: open(outwd+"out"+str(out)+"_"+x+".txt", "w"), pars)
for line in g:
  tmp = line.rstrip().split(" ")
  if len(tmp) < 2 :
    out += 1
    map(lambda x: x.close(), w0)
    w0 = map(lambda x: open(outwd+"out"+str(out)+"_"+x+".txt", "w"), pars)
  else :
    w0[pars.index(tmp[0])].write(" ".join(tmp[1:]) + "\n")

map(lambda x: x.close(), w0)
map(lambda x: os.remove(outwd+"out"+str(out)+"_"+x+".txt"), pars) # remove extra files




