sc = scatter(10*rand(100))

x = range(0,1,100)
y = x.^2

plot!(x,y,inset=(1,bbox(0.5,0.5,0.3,0.3)),subplot=2,xticks=0:1:1)