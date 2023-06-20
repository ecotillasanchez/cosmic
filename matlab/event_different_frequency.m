function [event,k,step,event1,event2,event3] = event_different_frequency(C,load_attack,t_max,frac_A,freq1,freq2,freq3)
% custom event setup for multiple frequency load oscillating attacks, as described in 
% F. Alanazi, J. Kim and E. Cotilla-Sanchez, "Load Oscillating Attacks 
% of Smart Grids: Vulnerability Analysis," in IEEE Access, vol. 11, 
% pp. 36538-36549, 2023, doi: 10.1109/ACCESS.2023.3266249.

% time_step=[1,20,10,5,(1/0.3),(1/0.6),(1/0.4),(1/0.5),(1/0.8),(1/0.7)];
time_step=[1,20,10,5,(1/0.3),(1/0.6),(1/0.4),(1/0.5),(1/0.8),(1/0.75),(1/0.9),(1/0.65),1/(2*0.42),1/(2*0.39)];
freq_all = randi(size(time_step),1,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EVENT MATRIX For AREA 1   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_attack_area_1 = load_attack(find (load_attack<=17));
if size(load_attack_area_1,2)==0
    load_attack_area_1 = unique(randi([1 17],1,1));
end
% freq = [7 5 2 3 1 2 3 4 1 7 2 8 4 3 2 3 5 3 5 3 7 2 4 5 4 8 3 6 5 4 5 4 2 6 7 1 2 7 6 2 3 4 2 5 4];
% freq = [7 3 2 2 2 
% freq = freq(itr); 
% freq = freq_all(1);
freq = freq1;
% load_frac=[0.1:0.1:0.9];
load_frac=frac_A(1,1);
frac = 1;
lembda_list=[0.1:0.1:1.2];

n_event = (((60/time_step(freq)))*length(load_attack_area_1))+2;
event_number=(60/time_step(freq));
if mod(60/time_step(freq),2) ~= 0
   n_event = (((60/time_step(freq))+1)*length(load_attack_area_1))+2;
   event_number=(60/time_step(freq))+1;
end
n_event = ceil(n_event);
event_1 = zeros(n_event,C.ev.cols);
% start
event_1(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip branches
k=1;
time=0;

for k=1:event_number
    if mod(k,2)~=0
        j=1;
        for i = (k-1)*length(load_attack_area_1)+2 : k*length(load_attack_area_1)+1
%       event_1(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_1(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_1(i,C.ev.shunt_loc) = load_attack_area_1(j);
%             event_1(i,C.ev.quantity) = ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_1(i,C.ev.quantity) = load_frac(frac);
            event_1(i,C.ev.change_by) = 1;
            j=j+1;
        end
        j=1;
    else
        for i = (k-1)*length(load_attack_area_1)+2 : k*length(load_attack_area_1)+1
%     event_1(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_1(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_1(i,C.ev.shunt_loc) = load_attack_area_1(j);
%             event_1(i,C.ev.quantity) = -ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_1(i,C.ev.quantity) = -load_frac(frac);
            event_1(i,C.ev.change_by) = 1;
            j=j+1;
    
        end
    end
    time=time+time_step(freq);
end
event1 = event_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EVENT MATRIX For AREA 2   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_attack_area_2 = load_attack(find (load_attack>17 & load_attack<=34));
if size(load_attack_area_2,2)==0
    load_attack_area_2 = unique(randi([18 34],1,1));
end
% freq = [4 3 7 5 1 4 4 1 3 4 4 6 4 4 8 3 8 4 3 4 7 3 1 7 2 3 5 3 2 2 5 6 5 4 5 2 4 3 4 3 3 8 4 7 3];
% freq = freq(itr); 
% freq = freq_all(2);
freq = freq2;
% load_frac=[0.1:0.1:0.9];
load_frac=frac_A(2,1);
frac = 1;
lembda_list=[0.1:0.1:1.2];

n_event = (((60/time_step(freq)))*length(load_attack_area_2))+2;
event_number=(60/time_step(freq));
if mod(60/time_step(freq),2) ~= 0
   n_event = (((60/time_step(freq))+1)*length(load_attack_area_2))+2;
   event_number=(60/time_step(freq))+1;
end
event_2 = zeros(n_event,C.ev.cols);
% start
event_2(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip branches
k=1;
time=0;

for k=1:event_number
    if mod(k,2)~=0
        j=1;
        for i = (k-1)*length(load_attack_area_2)+2 : k*length(load_attack_area_2)+1
%       event_2(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_2(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_2(i,C.ev.shunt_loc) = load_attack_area_2(j);
%             event_2(i,C.ev.quantity) = ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_2(i,C.ev.quantity) = load_frac(frac);
            event_2(i,C.ev.change_by) = 1;
            j=j+1;
        end
        j=1;
    else
        for i = (k-1)*length(load_attack_area_2)+2 : k*length(load_attack_area_2)+1
%     event_2(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_2(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_2(i,C.ev.shunt_loc) = load_attack_area_2(j);
%             event_2(i,C.ev.quantity) = -ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_2(i,C.ev.quantity) = -load_frac(frac);
            event_2(i,C.ev.change_by) = 1;
            j=j+1;
    
        end
    end
    time=time+time_step(freq);
end
event2 = event_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EVENT MATRIX For AREA 3   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_attack_area_3 = load_attack(find (load_attack>34 & load_attack<=51));
if size(load_attack_area_3,2)==0
    load_attack_area_3 = unique(randi([35 51],1,1));
end
% freq = [5 5 2 1 2 2 4 4 6 8 4 1 5 3 1 1 1 4 3 8 3 4 4 4 7 1 1 4 5 1 6 4 2 1 3 2 3 1 2 4 1 7 1 3 3]; 
% freq = freq(itr);
% freq = freq_all(3);
freq = freq3;
% load_frac=[0.1:0.1:0.9];
load_frac=frac_A(3,1);
frac = 1;
lembda_list=[0.1:0.1:1.2];

n_event = (((60/time_step(freq)))*length(load_attack_area_3))+2;
event_number=(60/time_step(freq));
if mod(60/time_step(freq),2) ~= 0
   n_event = (((60/time_step(freq))+1)*length(load_attack_area_3))+2;
   event_number=(60/time_step(freq))+1;
end
n_event = ceil(n_event);
event_3 = zeros(n_event,C.ev.cols);
% start
event_3(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip branches
k=1;
time=0;

for k=1:event_number
    if mod(k,2)~=0
        j=1;
        for i = (k-1)*length(load_attack_area_3)+2 : k*length(load_attack_area_3)+1
%       event_3(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_3(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_3(i,C.ev.shunt_loc) = load_attack_area_3(j);
%             event_3(i,C.ev.quantity) = ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_3(i,C.ev.quantity) = load_frac(frac);
            event_3(i,C.ev.change_by) = 1;
            j=j+1;
        end
        j=1;
    else
        for i = (k-1)*length(load_attack_area_3)+2 : k*length(load_attack_area_3)+1
%     event_3(i+1,[C.ev.time C.ev.type]) = [10+jj C.ev.trip_branch];
            event_3(i,[C.ev.time C.ev.type]) = [10+time C.ev.shed_load];
            event_3(i,C.ev.shunt_loc) = load_attack_area_3(j);
%             event_3(i,C.ev.quantity) = -ps_int.shunt(load_attack(j),2)*load_frac(frac);
            event_3(i,C.ev.quantity) = -load_frac(frac);
            event_3(i,C.ev.change_by) = 1;
            j=j+1;
    
        end
    end
    time=time+time_step(freq);
end
event3 = event_3;

event_time_A1 = 10;
event_time_A2 = 10;
event_time_A3 = 10;
pointer_area_1 = find(event_1(:,1)==event_time_A1);  
pointer_area_2 = find(event_2(:,1)==event_time_A2);
pointer_area_3 = find(event_3(:,1)==event_time_A3);

n_event = size(event_1,1)+size(event_2,1)-2+size(event_3,1)-2;
event =  zeros(n_event,C.ev.cols);
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
event_10=size(event_1(find(event_1(:,1)==event_time_A1),:),1)+size(event_2(find(event_2(:,1)==event_time_A2),:),1)+size(event_3(find(event_3(:,1)==event_time_A3),:),1);

event(2:event_10+1,:)=[event_1(find(event_1(:,1)==event_time_A1),:);event_2(find(event_2(:,1)==event_time_A2),:);event_3(find(event_3(:,1)==event_time_A3),:)];
event_pointer = event_10+2;

pointer_area_1 = pointer_area_1(end)+1;  
pointer_area_2 = pointer_area_2(end)+1;
pointer_area_3 = pointer_area_3(end)+1;

flag = NaN(3,1);
% for i = 1:(n_event-event_10)
while event_pointer <= size(event(:,1),1)-1
    
    pointer_area_1 = pointer_area_1(end);  
    pointer_area_2 = pointer_area_2(end);
    pointer_area_3 = pointer_area_3(end);
    if size(find(flag == inf)) ~= 0
%         pointer(flag == inf)= inf;
        z =flag == inf;
        if z(1)==1
            event_1(pointer_area_1)=inf;
        end
        if z(2) ==1
            event_2(pointer_area_2)=inf;
        end
        if z(3) ==1
            event_3(pointer_area_3)=inf;
        end
    end
    pointer = [event_1(pointer_area_1) event_2(pointer_area_2) event_3(pointer_area_3)];
    
   
    j=find(min(pointer)==pointer);
    if size(j,2)~=1
        j=j(1);
    end
    if j==1
        A1=size(find(event_1(:,1)==event_1(pointer_area_1,1)),1);
        pointer_area_1 = find(event_1(:,1)==event_1(pointer_area_1,1));
        d1=pointer_area_1(1);
        for d = event_pointer:event_pointer+A1-1
            event(d,:)= event_1(d1,:);
            d1=d1+1;
        end
        event_pointer=event_pointer+A1;
        if event_1(max(pointer_area_1)+1,1) ~=0
            pointer_area_1 = find(event_1(:,1)==event_1(max(pointer_area_1)+1,1));
        else
            flag(1,1) = inf;
        end
    elseif j==2
        A2=size(find(event_2(:,1)==event_2(pointer_area_2,1)),1);
        pointer_area_2 = find(event_2(:,1)==event_2(pointer_area_2,1));
        d2=pointer_area_2(1);
        for d = event_pointer:event_pointer+A2-1
            event(d,:)= event_2(d2,:);
            d2=d2+1;
        end
        event_pointer=event_pointer+A2;
        if event_2(max(pointer_area_2)+1,1) ~= 0
            pointer_area_2 = find(event_2(:,1)==event_2(max(pointer_area_2)+1,1));
        else
            flag(2,1) = inf;
        end
    elseif j ==3
        A3=size(find(event_3(:,1)==event_3(pointer_area_3,1)),1);
        pointer_area_3 = find(event_3(:,1)==event_3(pointer_area_3,1));
        d3=pointer_area_3(1);
        for d = event_pointer:event_pointer+A3-1
            event(d,:)= event_3(d3,:);
            d3=d3+1;
        end
        event_pointer=event_pointer+A3;
        if event_3(max(pointer_area_3)+1,1) ~= 0
            pointer_area_3 = find(event_3(:,1)==event_3(max(pointer_area_3)+1,1));
        else
            flag(3,1) = inf;
        end
        
    end
end
    
%     event(find(event_1(:,1)==10),:) = event_1(find(event_1(:,1)==10),:);



step =   time_step(freq);



event(end,[C.ev.time C.ev.type]) = [t_max C.ev.finish];
end