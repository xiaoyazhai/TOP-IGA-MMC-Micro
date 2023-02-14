function [tmpPhi]=tPhi(ComponentVariable,gridx,gridy,p)
st=ComponentVariable(7);
ct=sqrt(abs(1-st*st));
x1=ct*(gridx - ComponentVariable(1))+st*(gridy - ComponentVariable(2));
y1=-st*(gridx - ComponentVariable(1))+ct*(gridy - ComponentVariable(2));
bb=(ComponentVariable(5)+ComponentVariable(4)-2*ComponentVariable(6))/2/ComponentVariable(3)^2*x1.^2+(ComponentVariable(5)-ComponentVariable(4))/2*x1/ComponentVariable(3)+ComponentVariable(6);
tmpPhi= -((x1).^p/ComponentVariable(3)^p+(y1).^p./bb.^p-1);
end