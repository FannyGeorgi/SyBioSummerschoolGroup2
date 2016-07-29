clc; clear all; close all;
%parameter
tot=10;

%initialization
dboxes1=zeros(1,tot);
dboxes2=zeros(1,tot);
boxes=randi(10,1,tot);

%evaluating differences to the right
for i=1:tot
    if (i<tot)
        dboxes1(i)=boxes(i+1)-boxes(i);
    else
        dboxes1(i)=boxes(1)-boxes(i);
    end
end

%evaluating differences to the left
for i=1:tot
     if(i>1)
         dboxes2(i)=dboxes1(i-1);
     else
         dboxes2(1)=dboxes1(tot);
     end
end

disp('particles')
disp(boxes)
%disp('total number')
%disp(sum(boxes))  %check: mass conservation
disp('differnces')
disp(dboxes1)
disp(dboxes2)

%particle exchangers
pexchanger1=gradino(dboxes1);
pexchanger2=gradino(dboxes2);
%summing
%for i=1:tot; %right summing
%    if (i<tot)
%        boxes(i)=boxes(i)+dboxes(i)+dboxes(i-1);
%        boxes(i-1)=boxes(i-1)-dboxes(i);
%    else
%        boxes(i)=boxes(i)-dboxes(1)+dboxes(i-1);
%    end
%    disp(i)
%end

disp('final numbers')
disp(boxes)
disp(sum(boxes))  %check: mass conservation
