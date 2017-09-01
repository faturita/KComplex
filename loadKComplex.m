%function loadKComplexDataset()


[hdr, record] = edfread('C:\Users\User\Desktop\Data\KComplexes\excerpt1.edf');
Fs=200;

plot(record(3,round(50.4113*200):round(50.4113*200+0.6544*200)))

signal = record(3,:);
sg = bandpasseeg(signal', 1:1, Fs);

signal = sg;

%end


fid = fopen('C:\Users\User\Desktop\Data\KComplexes\Visual_scoring1_excerpt1.txt','r')

s=fscanf(fid, '%s', [1]);

events = fscanf(fid,'%f %f',[2 Inf]);
events=events'

fclose(fid)


for i=1:size(events,1)
    %figure;
    %plot(record(3,round(events(i,1)*Fs):round(events(i,1)*Fs+events(i,2)*Fs)))
    %figure;
    %plot(signal(1,round(events(i,1)*Fs)-200:round(events(i,1)*Fs+events(i,2)*Fs)+200))
    EEG(1,1,i).EEG = signal(round(events(i,1)*Fs):round(events(i,1)*Fs+1.2604*Fs),1);
    
end


imagescale=1;
subject=1;
minimagesize=250;
siftscale=[12 10];
siftdescriptordensity=1;
channelRange=1:1;

labelRange=ones(1,size(events,1));

epochRange=1:size(events,1);
for epoch=epochRange
    for channel=channelRange
        [eegimg, DOTS, zerolevel] = eegimage(channel,EEG(1,1,epoch).EEG,imagescale,1, false,minimagesize);
        label=labelRange(epoch);
        saveeegimage(subject,epoch,label,channel,eegimg);
        zerolevel = size(eegimg,1)/2;
        
        %             if ((size(find(trainingRange==epoch),2)==0))
        %                qKS=ceil(0.20*(Fs)*timescale):floor(0.20*(Fs)*timescale+(Fs)*timescale/4-1);
        %             else
        qKS=floor(size(eegimg,2)/2);
        %             end
        
        [frames, desc] = PlaceDescriptorsByImage(eegimg, DOTS,siftscale, siftdescriptordensity,qKS,zerolevel,false,'euclidean');
        
        
        F(channel,label,epoch).descriptors = desc;
        F(channel,label,epoch).frames = frames;
    end
end


for epoch=epochRange
    for channel=channelRange
        figure;DisplayDescriptorImageFull(F,1,epoch,1,1,-1,false);
    end
end

channel=1;
fprintf('Channel %d\n', channel);
fprintf('Building Test Matrix M for Channel %d:', channel);
[M, IX] = BuildDescriptorMatrix(F,channel,labelRange,epochRange);
fprintf('%d\n', size(M,2));


kdtree = vl_kdtreebuild(M) ;
fdist = [];

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

    K = size(M,2);
    distancetype='euclidean';
    
    [IDX, D] = vl_kdtreequery(kdtree,M,desc,'NumNeighbors',7);
    
    sumsrow = sum(sum(D));
    
    if (mod(slidewindow,200)==0)
        fprintf('.');
    end
    
    if (mod(slidewindow,200*10)==0)
        fprintf('\n%d:',slidewindow);
    end
            
    
%     [Z,I] = pdist2(M',desc',distancetype,'Smallest',K );
%     
%     k = 1;
% 
%     %Wi = Wdi(I(1:k,1:6)) ./ repmat( sum(Wdi(I(1:k,1:6))),k,1);     
%     Wi = ones(k,1);
%     if (k==1)
%         sumsrow = Z(1:k,1:1).*Wi(1:k,1:1);
%     else
%         sumsrow = dot(Z(1:k,1:1),Wi(1:k,1:1));
%     end
    
    fdist(end+1) = sumsrow;
    
end

