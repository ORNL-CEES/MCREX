a=poisson(2);
a=full(a);
load('z')
az=spconvert(z)
load('w')
aw=spconvert(w)
norm(eye(4)-az*a*aw')

