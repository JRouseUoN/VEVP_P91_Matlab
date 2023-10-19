time_split=[9007.042;
    18017.61;
    27031.69;
    36049.3;
    45070.42;
    54095.07;
    63123.24;
    72154.93;
    81190.14;
    90228.87;
    99271.13;
    108316.9];

data_proc=cell(32,1);

split_index=zeros(size(time_split,1),1);

for ii=1:1:size(time_split,1)
    split_index(ii,1)=find(abs(data(1:1083170,1)-time_split(ii))==min(abs((data(1:1083170,1)-time_split(ii)))));
end

for ii=1:1:size(time_split,1)
    if ii==1;
       data_proc{ii,1}=data(1:split_index(ii),:);
    else
        data_proc{ii,1}=data(split_index(ii-1):split_index(ii),:);
    end
end

data_temp=data(1083171:end,:);

data_temp_cyc=unique(data_temp(:,4));

jj=13;

for ii=1:1:size(data_temp_cyc,1)
    data_proc{jj,1}=data_temp(find(data_temp(:,4)==data_temp_cyc(ii)),:);
    jj=jj+1;
end
    

plotcolsvar=varycolor(size(data_proc,1));

hold on
for ii=1:1:size(data_proc,1)
    data_temp=data_proc{ii,1};
    data_temp(:,1)=data_temp(:,1)-data_temp(1,1);
    plot(data_temp(:,2),data_temp(:,3),'bx','MarkerEdgeColor',plotcolsvar(ii,:),'MarkerFaceColor',plotcolsvar(ii,:))
    clear data_temp
    pause
end
hold off
    
    
 
hold on
for ii=1:1:size(data_proc,1)
    data_temp=data_proc{ii,1};
    plot(data_temp(:,1),data_temp(:,3),'bx','MarkerEdgeColor',plotcolsvar(ii,:),'MarkerFaceColor',plotcolsvar(ii,:))
    clear data_temp
    pause
end
hold off











    