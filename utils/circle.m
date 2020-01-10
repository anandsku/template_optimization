function h = circle(x,y,r,ax,col,linestyle,linewidth)

th = 0:pi/50:2*pi;

xunit = r * cos(th) + x;

yunit = r * sin(th) + y;

h = plot(ax,xunit, yunit,'color',col,'linestyle',linestyle,'linewidth',linewidth);

