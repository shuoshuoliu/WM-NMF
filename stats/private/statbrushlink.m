function statbrushlink(handles,tf)
%STATBRUSHLINK Set brushing and linking for graphics objects
%   STATBRUSHLINK(H,TF) turns off brushing and linking for the objects
%   whose handles are given in the array H if TF is false, and turns it on
%   if TF is true.


%   Copyright 2007-2010 The MathWorks, Inc.

for j=1:numel(handles)
    gObj = handles(j);
    if ishghandle(gObj)
        brushBehavior = hggetbehavior(gObj,'Brush');
        brushBehavior.Enable = tf;
        linkBehavior = hggetbehavior(gObj,'Linked');
        linkBehavior.Enable = tf;
    end
end

