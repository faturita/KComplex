%function loadKComplexDataset()


[hdr, record] = edfread('C:\Users\User\Desktop\Data\KComplexes\excerpt1.edf');
Fs=200;

plot(record(3,round(50.4113*200):round(50.4113*200+0.6544*200)))

signal = record(3,:);
sg = bandpasseeg(signal', 1:1, Fs);

amplification = 12;
sg = zscore(sg) * amplification;

signal = sg;

%end


fid = fopen('C:\Users\User\Desktop\Data\KComplexes\Visual_scoring2_excerpt1.txt','r')

s=fscanf(fid, '%s', [1]);

events = fscanf(fid,'%f %f',[2 Inf]);
events=events'

fclose(fid)


for i=1:size(events,1)
    %figure;
    %plot(record(3,round(events(i,1)*Fs):round(events(i,1)*Fs+events(i,2)*Fs)))
    %figure;
    %plot(signal(1,round(events(i,1)*Fs)-200:round(events(i,1)*Fs+events(i,2)*Fs)+200))
    EEG(1,1,i).EEG = signal(max(round(events(i,1)*Fs-1.2600*Fs),1):round(events(i,1)*Fs+1.2600*2*Fs),1);
    
end


imagescale=1;
subject=1;
minimagesize=250;
siftscale=[12 10];
siftdescriptordensity=1;
channelRange=1:1;

labelRange=ones(1,size(events,1));

M = [];

epochRange=1:size(events,1);
for epoch=epochRange
    for channel=channelRange
        [eegimg, DOTS, zerolevel] = eegimage(channel,EEG(1,1,epoch).EEG,imagescale,1, false,minimagesize);
        label=labelRange(epoch);
        saveeegimage(subject,epoch,label,channel,eegimg);
        zerolevel = size(eegimg,1)/2;

        if (size(M,2)==0)
            qKS=floor(size(eegimg,2)/3);
        else
            qKS=floor(size(eegimg,2)/2)-200:floor(size(eegimg,2)/2)+200-1;
        end
        
        [frames, desc] = PlaceDescriptorsByImage(eegimg, DOTS,siftscale, siftdescriptordensity,qKS,zerolevel,false,'euclidean');
        
        
        F(channel,label,epoch).descriptors = desc;
        M = [M desc];
        F(channel,label,epoch).frames = frames;
    end
end

%%
for epoch=2:size(epochRange,2)
    DM = F(channel,label,epoch).descriptors;
    
    [ids, spdist] = knnsearch(DM',M(:,1)','k',1);
    
    F(channel,label,epoch).descriptors = F(channel,label,epoch).descriptors(:,ids);
    F(channel,label,epoch).frames = F(channel,label,epoch).frames(:,ids);
end
%%

for epoch=epochRange
    for channel=channelRange
        figure;DisplayDescriptorImageFull(F,1,epoch,1,1,-1,false);
    end
end


fesfsd

channel=1;
fprintf('Channel %d\n', channel);
fprintf('Building Test Matrix M for Channel %d:', channel);
[M, IX] = BuildDescriptorMatrix(F,channel,labelRange,epochRange);
fprintf('%d\n', size(M,2));


kdtree = vl_kdtreebuild(M) ;
fdist = [];

TM = [];

%%
for slidewindow=1:size(signal,1)-round(1.2604*Fs)-1
    signalsegment = signal(slidewindow:slidewindow+round(1.2604*Fs),1);    

    [eegimg, DOTS, zerolevel] = eegimage(channel,signalsegment,imagescale,1, false,minimagesize);

    %saveeegimage(subject,epoch,label,channel,eegimg);
    zerolevel = size(eegimg,1)/2;

    %             if ((size(find(trainingRange==epoch),2)==0))
    %                qKS=ceil(0.20*(Fs)*timescale):floor(0.20*(Fs)*timescale+(Fs)*timescale/4-1);
    %             else
    qKS=floor(size(eegimg,2)/2);
    %             end

    [frames, desc] = PlaceDescriptorsByImage(eegimg, DOTS,siftscale, siftdescriptordensity,qKS,zerolevel,false,'euclidean');
    TM = [TM desc];
    if (mod(slidewindow,200)==0)
        fprintf('.');
    end
    
    if (mod(slidewindow,200*10)==0)
        fprintf('\n%d:',slidewindow);
    end    
end


save('kcomplex2.mat');


%%
figure('Name','KComplex)','NumberTitle','off');
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
fcounter=1;
for i=1:size(M,2)
    ah=subplot_tight(6,5,fcounter,[0 0]);
    DisplayDescriptorImageFull(F,1,IX(i,3),IX(i,2),IX(i,1),IX(i,4),true);
    fcounter=fcounter+1;
end

%%
% for i=1:size(TM,2)
%     K = size(M,2);
%     distancetype='euclidean';
%     
%     [IDX, D] = vl_kdtreequery(kdtree,M,desc,'NumNeighbors',1);
%     
%     sumsrow = sum(sum(D));
%     
%     if (mod(slidewindow,200)==0)
%         fprintf('.');
%     end
%     
%     if (mod(slidewindow,200*10)==0)
%         fprintf('\n%d:',slidewindow);
%     end
%             
%     
% %     [Z,I] = pdist2(M',desc',distancetype,'Smallest',K );
% %     
% %     k = 1;
% % 
% %     %Wi = Wdi(I(1:k,1:6)) ./ repmat( sum(Wdi(I(1:k,1:6))),k,1);     
% %     Wi = ones(k,1);
% %     if (k==1)
% %         sumsrow = Z(1:k,1:1).*Wi(1:k,1:1);
% %     else
% %         sumsrow = dot(Z(1:k,1:1),Wi(1:k,1:1));
% %     end
%     
%     fdist(end+1) = sumsrow;
%     
% end

%%
kdtree = vl_kdtreebuild(M(:,1:15)) ;

TTM = TM(:,round(1.26*Fs):end);

[IDX, D] = vl_kdtreequery(kdtree,M(:,1:15),TM,'NumNeighbors',7);

if (size(D,1) > 1)
    % Sumarizar los valores en la primera dimension de manera de
    D = sum(D,1);
end

fdist=D;
[values,order] = sort(fdist);

eventssample = round(events(:,1)*Fs+1.2604/2*Fs);


tp=zeros(1,size(eventssample,1));
fp=[];
matchingdescriptorsids=[];
nomatchingdescriptorsids=[];

S = size(eventssample,1);

S=size(order,2);

for i=1:S
    
    [ids, spdist] = knnsearch(eventssample,order(i),'k',1);
    
    if (spdist<1.26*Fs)
        tp(ids) = 1;
        matchingdescriptorsids(end+1) = i;
    else
        nomatchingdescriptorsids(end+1) = i;
    end
    
    if (size(find(tp==0),2)==0)
        break;
    end
end

i
tp

figure;
hold on;
for a=1:size(eventssample,1)
    plot(eventssample(a),1,'rx');
end
for a=1:size(matchingdescriptorsids,2)
   plot(order( matchingdescriptorsids(a)),2,'b.');
end
for a=1:size(nomatchingdescriptorsids,2)
   plot(order( nomatchingdescriptorsids(a)),2,'r.');
end
axis([0 360000 0 10])
hold off;

orderevents=[];

for a=1:i
    orderevents(end+1) = order(a);
end

%%
predicted=zeros(1,size(signal,1));
expected=zeros(1,size(signal,1));

predicted(orderevents)=1;
expected(eventssample)=1;

pred=[];
expe=[];

for slidewindow=1:Fs:size(signal,1)-round(1.2604*Fs)-1
    if ((size(find(predicted(slidewindow:slidewindow+Fs)==1),2))>0)
        pred(end+1)=1;
    else
        pred(end+1)=0;
    end
    
    if ((size(find(expected(slidewindow:slidewindow+Fs)==1),2))>0)
        expe(end+1)=1;
    else
        expe(end+1)=0;
    end    
        
end

C=confusionmat(expe,pred)

ACC = (C(1,1)+C(2,2)) / size(pred,2)

tpr = (C(2,2) / (C(2,2)+C(2,1)))
fdr = (C(1,2) / (C(1,2)+C(1,1)))


%%

figure;
hold on;
plot(fdist);
for i=1:size(events,1)
    %figure;
    %plot(record(3,round(events(i,1)*Fs):round(events(i,1)*Fs+events(i,2)*Fs)))
    %figure;
    %plot(signal(1,round(events(i,1)*Fs)-200:round(events(i,1)*Fs+events(i,2)*Fs)+200))
    %EEG(1,1,i).EEG = signal(round(events(i,1)*Fs):round(events(i,1)*Fs+1.2604*Fs),1);
    
    plot( round(events(i,1)*Fs+1.2604/2*Fs),0.5*10^6,'rx');
    plot( order(i),0.5*10^6,'b.');
end
hold off;


