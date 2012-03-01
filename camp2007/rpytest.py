import rpy,random

x = [random.random() for i in range(100)]
y = [random.random() for i in range(100)]

s = 'short'
xlab='x'
ylab='y'
print 'title=',s
rpy.r.pdf('short.pdf',h=3,w=3)
rpy.r.plot(x,y,pch='1',col='red',xlab=xlab,ylab=ylab,main=s)
rpy.r.dev_off()
rpy.r.pdf('long.pdf',h=3,w=3)
s= 'a very, very, very long main title and getting even longer'
print 'title=',s
rpy.r.plot(x,y,pch='1',col='red',xlab=xlab,ylab=ylab,main=s)
rpy.r.dev_off()

