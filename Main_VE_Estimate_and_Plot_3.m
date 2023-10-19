

hold on
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{1,1};
    if ii==1
        plot(data_temp(:,1),data_temp(:,3),'rx')
    elseif ii==2
        plot(data_temp(:,1),data_temp(:,3),'bx')
    elseif ii==3
        plot(data_temp(:,1),data_temp(:,3),'gx')
    end
    clear data_temp
end
hold off