function y = correctLabel(yz, fd)
yz = castLabel(yz,-1);
yz(fd==1) = yz(fd==1) * -1;
y = yz;