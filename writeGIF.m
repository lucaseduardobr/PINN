function writeGIF(hax,gname,opt)

if nargin == 2
    opt = 1;
end

if opt == 0
  frame = getframe(hax); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  imwrite(imind,cm,gname,'gif', 'Loopcount',inf); 
else
  frame = getframe(hax); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  imwrite(imind,cm,gname,'gif','WriteMode','append');     
end