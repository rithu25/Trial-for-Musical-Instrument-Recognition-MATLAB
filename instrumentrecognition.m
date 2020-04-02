function instrumentrecognition(val)
music=strcat('C:\Users\prashaanth\Desktop\audio\FASST_version_v1\test\',val,'.wav');
J_spat=instrumentclassification(music);
file={'C:\Users\prashaanth\Desktop\audio\FASST_version_v1\output_data\'};
load 'mfcccombininghor.mat';
instrumentscell={'cello','clarinet','elecguitar','flute','gac','organ','piano','sax','trumpet','tuba','violin'};
mfcccell=cell(1,J_spat);
% for i=1:11
%     mfccmatcell=cell(1,10);
%     temp2=strcat(directory,instrumentscell{i});
%     for n = 1:10
%         fprintf('%d of %d\n', n, 10);
%         if(i~=7)
%             musicfiles=strcat(temp2,'\','filename_ (',num2str(n),').aif');
%             
%         else
%             musicfiles=strcat(temp2,'\','filename_ (',num2str(n),').aiff');
%         end
%         a=mirframe(musicfiles,0.25);
%         mfccmatcell{n}=mirgetdata(mirmfcc(a));
%     end
%     mfcctotalcell{i}=cell2mat(mfccmatcell);
% end
% mfcccell=cell(1,11);
% for i=1:11
%     mfcctotal=mfcctotalcell{i};
%     [rows,columns]=size(mfcctotal);
%     mfccmat=zeros(rows,1);
%     for i1=1:rows
%         value=0;
%         for i2=1:columns
%             value=value+mfcctotal(i1,i2);
%         end
%         mfccmat(i1)=value/columns;
%     end
%     mfcccell{i}=mfccmat;
% end
% for j=1:11
%     mfcccombiningcell{j}=cat(2,mfcctotalcell{1,j}{:});
% end
for i=1:J_spat
    temp1=mirframe(strcat(file,'output',num2str(i),'.wav'),0.25);
    mfcc=mirgetdata(mirmfcc(temp1));
    mfcccell{i}=mfcc;
end
target=zeros(11,10);
for i=1:10
target(:,i)=[1 zeros(1,10)];
end 
for i=11:20
target(:,i)=[zeros(1,1) 1 zeros(1,9)];
end
for i=21:30
target(:,i)=[zeros(1,2) 1 zeros(1,8)];
end
for i=31:40
target(:,i)=[zeros(1,3) 1 zeros(1,7)];
end
for i=41:50
target(:,i)=[zeros(1,4) 1 zeros(1,6)];
end
for i=51:60
target(:,i)=[zeros(1,5) 1 zeros(1,5)];
end
for i=61:70
target(:,i)=[zeros(1,6) 1 zeros(1,4)];
end
for i=71:80
target(:,i)=[zeros(1,7) 1 zeros(1,3)];
end
for i=81:90
target(:,i)=[zeros(1,8) 1 zeros(1,2)];
end
for i=91:100
target(:,i)=[zeros(1,9) 1 zeros(1,1)];
end
for i=101:110
target(:,i)=[zeros(1,9) 1 zeros(1,1)];
end
% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainlm';  % Scaled conjugate gradient backpropagation.

% Create a Pattern Recognition Network
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize, trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
% [net,~] = train(net,mfcccombininghor,target);

% Test the Network
% y = net(mfcccombininghor);
% View the Network
% view(net)
finaltext=cell(1,J_spat);
finaltextcell=cell(1,J_spat);
for i=1:J_spat
    [Y,~,~]=myNeuralNetworkFunction1(mfcccell{i});
    % Plots
    % Uncomment these lines to enable various plots.
    %figure, plotperform(tr)
    %figure, plottrainstate(tr)
    %figure, ploterrhist(e)
    %figure, plotconfusion(t,y)
    %figure, plotroc(t,y)
    [row,col]=size(Y);
    indexmat=zeros(1,col);
    for p=1:col
        for q=1:row
            ymat=Y(:,p);
            ymax=max(ymat);
            if(isequal(ymax,Y(q,p)))
                break;
            end
        end
        indexmat(p)=q;
    end
    indexunique=unique(indexmat);
    [~,cols]=size(indexunique);
    cntmat=zeros(1,cols);
    for r=1:cols
        cnt=0;
        for s=1:col
            if(isequal(indexunique(r),indexmat(s)))
                cnt=cnt+1;
            end
        end
        cntmat(r)=cnt;
    end
    cntmax=max(cntmat);
    finalindex=zeros(1,cols);
    for t=1:cols
        if(isequal(cntmax,cntmat(1,t)))
            finalindex(t)=indexunique(t);
        end
    end
    textcell=cell(1,cols);
    for r=1:cols
        for s=1:11
            text='';
            if(isequal(finalindex(r),s))
                text=strcat(instrumentscell{s},' is played');
                break;
            end
        end
        textcell{r}=text;
    end
    finaltext{i}=textcell;
end
for i=1:J_spat
    final=unique(finaltext{i});
    [~,cols]=size(final);
    for j=1:cols
        if(~isempty(final{j}))
            finaltextcell{i}=final{j};
            break;
        end
    end
end
uniquefinaltext=unique(finaltextcell);
[~,col]=size(uniquefinaltext);
for i=1:col
    msgbox(uniquefinaltext{i});
end