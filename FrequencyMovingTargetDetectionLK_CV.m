clc;
clear;
%addpath(genpath(pwd));
addpath(genpath('myfun'));
%addpath('myfun');
dbstop if error;

TimeHeight=11;%19;
deltay=zeros(TimeHeight-1,1);
deltarow=zeros(TimeHeight,1);
TimeCenter=(TimeHeight+1)/2;
TimeCenter=floor(TimeCenter);
%mov=mmreader('D:\视频资料\PTZCar320.avi'); %读入视频
%mov=mmreader('D:\视频资料\ParkingLot2Seg320.avi'); %Road1Seg ParkingLot2Seg320 MyVideo320
%fnum=mov.NumberOfFrames; %读取视频的祯数，mov为1*temp
%FrameHeight=mov.height;%240;
%FrameWidth=mov.width;%320;

%FirstFrame=50; fnum=250; dataName='conference1';

% FirstFrame=15; fnum=180; dataName='ParkingLot1Seg';
% FirstFrame=15; fnum=180; dataName='ParkingLot2Seg';
% FirstFrame=75; fnum=260; dataName='MyVideo320';
% FirstFrame=20; fnum=230; dataName='Road1Seg';

% FirstFrame=740; fnum=745; dataName='continuousPan';

%D:\视频资料\PTZ\continuousPan\input
% filenamestr=['D:\视频资料\PTZ\continuousPan\input\in'];
% filename=sprintf('%s%06d.jpg',filenamestr,FirstFrame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mov=mmreader('D:\视频资料\Lab\5.AVI');
% fnum=mov.NumberOfFrames;
% for f=1:fnum
%     A=read(mov,f);
%     A=imresize(A,[240 320]);
%     filename=sprintf('%s%d.bmp','D:\视频资料\Lab\5\',f);
%     imwrite(A,filename);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataName='Lab1_1'; %'MyVideo320'; %'Lab4_1'; %'ParkingLot1Seg';%'ParkingLot2Seg';%'Road1Seg'; %

switch dataName
    case 'Road1Seg'
        FirstFrame=20; fnum=230; totalframe=210;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\RoadTrueRect.txt';
    case 'MyVideo320'
        FirstFrame=75; fnum=260; totalframe=185;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\ConferenceTrueRect.txt';
    case 'ParkingLot1Seg'
        FirstFrame=15; fnum=180; totalframe=165;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\ParkingLot1TrueRect.txt';
    case 'ParkingLot2Seg'
        FirstFrame=15; fnum=180; totalframe=165;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\ParkingLot2TrueRect.txt';
    case 'Lab1_1'
        FirstFrame=290; fnum=370; totalframe=80;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\Lab1_1TrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\1\Seg1_RightDown\'];
    case 'Lab4_1'
        FirstFrame=160; fnum=310; totalframe=150;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\Lab4_1TrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\4\Seg1_Down\'];
        
    case 'Lab1_2'
        FirstFrame=2010; fnum=2160;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\1\Seg2_LeftTop\'];
    case 'Lab1_3'
        FirstFrame=2680; fnum=2800;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\1\Seg3_RightDown\'];
    case 'Lab2_1'
        FirstFrame=325; fnum=475;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\2\Seg1_Left\'];
    case 'Lab2_2'
        FirstFrame=910; fnum=1010;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\2\Seg2_Right\'];
    case 'Lab5_1'
        FirstFrame=160; fnum=350;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\5\Seg1_LeftTop\'];
    case 'Lab5_2'
        FirstFrame=820; fnum=980;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\LabTrueRect.txt';
        %filenamestr1=['D:\视频资料\Lab\5\Seg2_RightDown\'];

    case 'OSU_Color_2b'
        FirstFrame=15; fnum=371;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\ParkingLot2TrueRect.txt';
    case 'Walk1'
        %FirstFrame=40; fnum=170;
        FirstFrame=270; fnum=470;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\ParkingLot2TrueRect.txt';
    otherwise
        FirstFrame=20; fnum=230;
        fid1='E:\TeacherMonocular\TLD_source\resultimg\RoadTrueRect.txt';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenamestr=['D:\视频资料\' dataName '\'];
% for f=FirstFrame:fnum
%     filename=sprintf('%s%d.bmp',filenamestr1,f);
%     A=imread(filename);
%     filename=sprintf('%s%05d.bmp',filenamestr,f);
%     imwrite(A,filename);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenamestr=['D:\视频资料\' dataName '\'];
filename=sprintf('%s%05d.bmp',filenamestr,FirstFrame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ppm转bmp
% filenamestr1=['F:\pami2013_MoSeg\Results_Total\Lab4_1\SparseSegmentation\Segments'];%1\marple8_'];
% FirstFrame=0; fnum=150;
% for f=FirstFrame:fnum
%     filename=sprintf('%s%03d.ppm',filenamestr1,f);
%     A=imread(filename);
%     filename=sprintf('%s%03d.bmp',filenamestr1,f);
%     imwrite(A,filename);
% end
% a=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%视频序列分段
% filenamestr1=['F:\pami2013_MoSeg\moseg\TrainingSet\Data\ourdata\' dataName '\'];%1\marple8_'];
% seg=fnum-FirstFrame+1;
% for f=FirstFrame:fnum %25:fnum %
%     filename=sprintf('%s%05d.bmp',filenamestr,f);
%     A=imread(filename);
%     frame=f-FirstFrame;
%     docname=floor(frame/seg)+1;
%     framename=mod(frame,seg)+1;
%     filename=sprintf('%s%d%s%02d.jpg',filenamestr1,docname,'\marple8_',framename);
%     imwrite(A,filename);
% end
% a=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%LongTermTtajectory Precision
% % fid1=['E:\TeacherMonocular\TLD_source\resultimg\ConferenceTrueRect.txt'];
% % fid1=['E:\TeacherMonocular\TLD_source\resultimg\ParkingLot1TrueRect.txt'];
% % fid1=['E:\TeacherMonocular\TLD_source\resultimg\ParkingLot2TrueRect.txt'];
% % fid1=['E:\TeacherMonocular\TLD_source\resultimg\RoadTrueRect.txt'];
% text1=fopen(fid1,'rt');
% fid2='resultimg\SegF1.txt';
% text2=fopen(fid2,'wt');
% filenamestr1=['F:\pami2013_MoSeg\Results_bmp\total\' dataName '\Segments'];
% %filenamestr1=('F:\pami2013_MoSeg\Results_bmp\total\Lab1_1\Segments');
% filename1=sprintf('%s%03d.bmp',filenamestr1,1);
% A=imread(filename1);
% A_Motion=zeros(size(A,1),size(A,2));
% A_Des=A_Motion;
%     
% %FirstFrame=0; fnum=80;
% for f=0:totalframe
%     filenamestr1=['F:\pami2013_MoSeg\Results_bmp\total\' dataName '\Segments'];
%     %filenamestr1=('F:\pami2013_MoSeg\Results_bmp\total\Lab1_1\Segments');
%     filename1=sprintf('%s%03d.bmp',filenamestr1,f);
%     A=imread(filename1);
%     A=ChangeImg(A);
%     [A_Motion A_Hist]=MainStat(A);
% 
%     LeftX=fscanf(text1,'%f',1);
%     TopY=fscanf(text1,'%f',1);
%     RightX=fscanf(text1,'%f',1);
%     BotY=fscanf(text1,'%f',1);
%     GT=(RightX-LeftX+1)*(BotY-TopY+1);
%     A_Patch=A(TopY:BotY,LeftX:RightX,:);
%     A_Patch=ChangeImg(A_Patch);
%     [A_Patch_Motion A_Patch_Hist]=MainStat(A_Patch);
%     maxP=-1; maxP_value=0; TP=0; ForeTotal=0;
%     for i=1:size(A_Patch_Hist,1)
%         value=A_Patch_Hist(i,1);
%         TP_tmp=A_Patch_Hist(i,2);
%         ForeTotal_tmp=A_Hist(A_Hist(:,1)==value,2);
%         P=TP_tmp/ForeTotal_tmp;
%         if P>maxP
%             TP=TP_tmp;
%             ForeTotal=ForeTotal_tmp;
%             maxP=P;
%             maxP_value=value;
%         end
%     end
%     tmp=find(A_Motion==maxP_value);
%     A_Des(:)=0;
%     A_Des(tmp)=255;
%     filenamestr1=['F:\pami2013_MoSeg\Results_bmp\total\' dataName '\SegBinary\'];
%     %filenamestr1=('F:\pami2013_MoSeg\Results_bmp\total\Lab1_1\SegBinary\');
%     filename1=sprintf('%s%d.bmp',filenamestr1,f);
%     imwrite(A_Des,filename1);
% 
%     filename=sprintf('%s%05d.bmp',filenamestr,FirstFrame+f);
%     A=imread(filename);
%     A1=A(:,:,1);
%     tmp=find(A_Des>0);
%     A1(tmp)=255;
%     A(:,:,1)=A1;
%     filenamestr1=['F:\pami2013_MoSeg\Results_bmp\total\' dataName '\SegCol\'];
%     %filenamestr1=('F:\pami2013_MoSeg\Results_bmp\total\Lab1_1\SegCol\');
%     filename1=sprintf('%s%d.bmp',filenamestr1,f);
%     imwrite(A,filename1);
%     
%     P=TP/ForeTotal;
%     R=TP/GT;
%     F1=2*TP/(GT+ForeTotal);
%     fprintf(text2,'%f ',TP);
%     fprintf(text2,'%f ',ForeTotal);
%     fprintf(text2,'%f ',GT);
%     fprintf(text2,'%f ',P);
%     fprintf(text2,'%f ',R);
%     fprintf(text2,'%f\n',F1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=imread(filename);
CellA=cell(TimeHeight,1);
[FrameHeight,FrameWidth,~]=size(A);
EdgeDesTmp=zeros(TimeHeight,FrameWidth);
C=zeros(FrameHeight,FrameWidth,'uint8');
CPre=zeros(FrameHeight,FrameWidth,'uint8');
CTotal=zeros(FrameHeight,FrameWidth,'uint8');
CTotalPre=255*ones(FrameHeight,FrameWidth,'uint8');
CPre2=zeros(FrameHeight,FrameWidth,'uint8');
C1=zeros(FrameHeight,FrameWidth,'uint8');
Target=zeros(FrameHeight,FrameWidth,3,'uint8');
APre=zeros(FrameHeight,FrameWidth,'uint8');
SelectFlag=0;
Quan=1;
slope=zeros(FrameHeight,Quan);
Error=0;
MinPos=0;
MinError=1000;
ObjH=0;ObjW=0;
% WholeROI.TopLeftX=20;WholeROI.BotRightX=300;
% WholeROI.TopLeftY=50;WholeROI.BotRightY=195;
WholeROI.TopLeftX=2;WholeROI.BotRightX=FrameWidth-1;
WholeROI.TopLeftY=2;WholeROI.BotRightY=FrameHeight-1;
ROI=WholeROI;
margin=20; %10; %
WholeTimeROI.TopLeftX=margin;WholeTimeROI.BotRightX=FrameWidth-margin;
WholeTimeROI.TopLeftY=1;WholeTimeROI.BotRightY=TimeHeight;
TimeROI=WholeTimeROI;
FrameIndex=1;
FrameIndex1=1;
se=strel('square',2);
warning off;
lk(0);
BB1=[2;2;FrameWidth-1;FrameHeight-1];
xFI=bb_points(BB1,100,100,20);%bb_points(BB1,10,10,5); % generate 10x10 grid of points within BB1 with margin 5 px
%FirstFrame=25;
for i=TimeHeight-1:-1:1
        filename=sprintf('%s%05d.bmp',filenamestr,FirstFrame-i);
        
%         filenamestr=['D:\视频资料\PTZ\continuousPan\input\in'];
%         filename=sprintf('%s%06d.jpg',filenamestr,FirstFrame-i);
        A=imread(filename);
        %A=read(mov,FirstFrame-i);
        A=rgb2gray(A);
        if i<TimeHeight-1
            A1=CellA{FrameIndex-1,1};
            xFJ    = lk(2,A,A1,xFI,xFI); % track all points by Lucas-Kanade tracker from frame I to frame J, estimate Forward-Backward error, and NCC for each point
            medFB  = median2(xFJ(3,:)); % get median of Forward-Backward error
            medNCC = median2(xFJ(4,:)); % get median for NCC
            idxF   = xFJ(3,:) <= medFB & xFJ(4,:)>= medNCC; % get indexes of reliable points
            [Dx Dy maxHistDy maxHistDyid delta1]=bb_predict(BB1,xFI(:,idxF),xFJ(1:2,idxF)); % estimate BB2 using the reliable points only
            Dx=round(Dx);Dy=round(Dy);
%             [Dx, Dy]=Lucas_Kanade(A,A1,2);%Img1=Img2+D
%             mindy=min(Dy(:));
%             maxdy=max(Dy(:));
%             lengthHistDy=maxdy-mindy+1;
%             HistDy=hist(Dy(:),lengthHistDy);
%             [maxHistDy maxHistDyid]=max(HistDy);
%             deltay(FrameIndex1)=maxHistDyid+mindy-1;
            deltay(FrameIndex1)=Dy;
            FrameIndex1=FrameIndex1+1;
            if FrameIndex1>TimeHeight-1
                FrameIndex1=1;
            end
        end
        CellA(FrameIndex,1)={A};
        FrameIndex=FrameIndex+1;
        if FrameIndex>TimeHeight
            FrameIndex=1;
        end
end
% fid1=['E:\TeacherMonocular\TLD_source\resultimg\RoadTrueRect.txt'];
% fid1=['E:\TeacherMonocular\TLD_source\resultimg\ConferenceTrueRect.txt'];
% fid1=['E:\TeacherMonocular\TLD_source\resultimg\ParkingLot1TrueRect.txt'];
% fid1=['E:\TeacherMonocular\TLD_source\resultimg\ParkingLot2TrueRect.txt'];
text1=fopen(fid1,'rt');
fid2='resultimg\SegF1.txt';
text2=fopen(fid2,'wt');
fid3='resultimg\Segtime.txt';
text3=fopen(fid3,'wt');
for f=FirstFrame:fnum %25:fnum %
    %tic
    B=zeros(TimeHeight,FrameWidth,'uint8');%行序列
    %A=read(mov,f);
    
    filename=sprintf('%s%05d.bmp',filenamestr,f);
    
%     filenamestr=['D:\视频资料\PTZ\continuousPan\input\in'];
%     filename=sprintf('%s%06d.jpg',filenamestr,f);

    A=imread(filename);
    %A=imread(filename);
    tic
    A=rgb2gray(A);
    %A=imresize(A,[FrameWidth FrameWidth]);
    CellA(FrameIndex,1)={A};
    FrameIndexTmp=FrameIndex-1;
    if FrameIndexTmp<1
       FrameIndexTmp=FrameIndexTmp+TimeHeight; 
    end
    A1=CellA{FrameIndexTmp,1};
    xFJ    = lk(2,A,A1,xFI,xFI); % track all points by Lucas-Kanade tracker from frame I to frame J, estimate Forward-Backward error, and NCC for each point
    medFB  = median2(xFJ(3,:)); % get median of Forward-Backward error
    medNCC = median2(xFJ(4,:)); % get median for NCC
    idxF   = xFJ(3,:) <= medFB & xFJ(4,:)>= medNCC; % get indexes of reliable points
    [Dx Dy maxHistDy maxHistDyid delta1]= bb_predict(BB1,xFI(:,idxF),xFJ(1:2,idxF)); % estimate BB2 using the reliable points only
    Dx=round(Dx);Dy=round(Dy);
%     [Dx, Dy]=Lucas_Kanade(A,A1,2);%Img1=Img2+D
%     mindy=min(Dy(:));
%     maxdy=max(Dy(:));
%     lengthHistDy=maxdy-mindy+1;
%     HistDy=hist(Dy(:),lengthHistDy);
%     [maxHistDy maxHistDyid]=max(HistDy);
%     deltay(FrameIndex1)=maxHistDyid+mindy-1;
    deltay(FrameIndex1)=Dy;
    
    if(~SelectFlag)
        CPre2(:)=0;
    elseif(SelectFlag==1)
        CPre(:)=0;
    else
        C(:)=0;
    end
    FrameIndexTmp1=FrameIndex1;
    deltarow(1)=0;
    for i=2:TimeHeight
       deltarow(i)=deltay(FrameIndexTmp1);
       FrameIndexTmp1=FrameIndexTmp1-1;
       if FrameIndexTmp1<1
           FrameIndexTmp1=FrameIndexTmp1+TimeHeight-1;
       end
    end
   
    for k=ROI.TopLeftY:ROI.BotRightY %75:ROI.BotRightY %104:ROI.BotRightY %k=130 %k=50:240 %k=110:240 %145 %75%108
        FrameIndexTmp=FrameIndex;
        row=k;
        for i=1:TimeHeight
            A=CellA{FrameIndexTmp,1};
            row=row-deltarow(i);
            row=max(WholeROI.TopLeftY,row);
            row=min(WholeROI.BotRightY,row);
            B(i,:)=A(row,:);
            FrameIndexTmp=FrameIndexTmp-1;
            if FrameIndexTmp<1
                FrameIndexTmp=FrameIndexTmp+TimeHeight;
            end
        end
    %figure(1);
    %imshow(B);
%       string='D:\Visual Tracking\Frequency moving target detection\DynamicBG MovTargetDetect\resultimg\';
%       string=strcat(string,[num2str(f),'PTZCar.bmp']);
%       imwrite(A,string);
%   string='D:\Visual Tracking\Frequency moving target detection\DynamicBG MovTargetDetect\resultimg\';
%       string=strcat(string,[num2str(f),'PTZCarTime.bmp']);
%       imwrite(B,string);

        Edge=edge(B,'canny',0.2);%edge(B,'canny',0.13); %
%         figure(1);
%         imshow(Edge);
%         imwrite(Edge,'F25_R104_Edge.bmp');
        Edge(:,1:TimeROI.TopLeftX)=0;
        Edge(:,TimeROI.BotRightX:FrameWidth)=0;
        FiltEdge=EdgeFilt(Edge,TimeROI,3);
%         figure(2);
%         imshow(FiltEdge);
%         imwrite(FiltEdge,'F25_R104_FiltEdge.bmp');
        [Slope,EdgeDes,EdgeDesNum,EdgeNum]=CalcEdgeSlope(FiltEdge,TimeROI);
        
%      string='D:\Visual Tracking\Frequency moving target detection\DynamicBG MovTargetDetect\resultimg\';
%       string=strcat(string,[num2str(f),'PTZCarEdgeFilt.bmp']);
%       imwrite(FiltEdge,string);

        Slope(isnan(Slope)) = [];
        
        if ~isempty(Slope)
            %[inliers,idx,outliers] = deleteoutliers(Slope);
            inliers=Slope;
            %inliers(inliers==20)=[];
            inliers(inliers==5|inliers==20)=[];
            if ~isempty(inliers)
                delta=abs(mean(abs(inliers))/4);%6;
            else
                delta=abs(mean(abs(Slope))/6);
            end
        end
        %Slope(abs(Slope)>=3&abs(Slope)<20)=5;
        Slope(abs(Slope)>=5&abs(Slope)<20)=5;
        %tmp=find(abs(Slope)>=3);
        %Slope(tmp)=5;
        if ~isempty(Slope)
            [y,center,classnum] = HierarchicalCluster(Slope,delta);
            %[y,center,classnum]=ClusterEdge(Slope,delta);
            g=y(:,2);
            radius=zeros(length(classnum),1);
            ColVar=radius;
            classnumTmp=classnum;
            classnumTmp(center==20)=0;
            [maxvalue maxclass]=max(classnumTmp);
            for i=1:length(classnum)
                if classnum(i)~=0
                    tmp=find(g==i);
                    if ~isempty(tmp)
                        if classnum(i)>0.6*maxvalue %0.4*maxvalue %
                            Slope(tmp)=20;
                        else
                            EdgeDesTmp(:)=0;
                            for j=1:length(tmp)
                                value=tmp(j)+EdgeNum-1;
                                tmp1=find(EdgeDes==value);
                                if length(tmp1)<6
                                    EdgeDes(tmp1)=0;
                                    EdgeDesNum(EdgeDesNum==value)=[];
                                else
                                    EdgeDesTmp(tmp1)=EdgeDes(tmp1);
                                end
                            end
                            [radius(i) ColVar(i)]=CalcRadius(EdgeDesTmp);
                        end
                    end
                end
            end
            %tmp=find((radius>originalwidth/2));%|(radius<originalwidth/10));
            tmp=find(ColVar>0.1);
            if ~isempty(tmp)
                for i=1:length(tmp)
                    tmp1=find(g==tmp(i));
                    Slope(tmp1)=20;
                end
            end
        end
        DesSlope=find(Slope==20);
        num=size(DesSlope);
        if ~isempty(DesSlope)
            for i=1:num(1)
                value=EdgeNum+DesSlope(i)-1;
                tmp= EdgeDes==value;
                EdgeDes(tmp)=0;
                EdgeDesNum(EdgeDesNum==value)=[];
            end
        end
%     figure(3);
%          imshow(EdgeDes);
         %imwrite(EdgeDes,'F25_R104_EdgeDes.bmp');
        if ~isempty(EdgeDesNum)
            DesCol=[];
            for i=1:length(EdgeDesNum)
                [Row Col]=find(EdgeDes==EdgeDesNum(i));
                [minRow,minRowId]=min(Row);
                index=EdgeDesNum(i)-EdgeNum+1;
                if Slope(index)>0 && Slope(index)~=20
                    DesCol=[DesCol Col(minRowId)-minRow-2:Col(minRowId)-minRow+4];
                elseif Slope(index)<0
                    DesCol=[DesCol Col(minRowId)+minRow-4:Col(minRowId)+minRow+2];
                else
                    DesCol=[DesCol Col(minRowId)-3:Col(minRowId)+3];
                end
            end
            %DesCol=[DesCol-3:DesCol+3];
            %DesCol=[DesCol;DesCol-1;DesCol-2;DesCol+1;DesCol+2];
            DesCol=max(DesCol,ROI.TopLeftX);
            DesCol=min(DesCol,ROI.BotRightX);
            
            if(~SelectFlag)
                CPre2(k,DesCol)=255;
            elseif(SelectFlag==1)
                CPre(k,DesCol)=255;
            else
                C(k,DesCol)=255;
            end
        end
    end
    
    if(~SelectFlag)
        CPre2=FiltImg(CPre2,6);
        CPre2=FiltImg(CPre2,4);
        CPre2=imdilate(CPre2,se);
    elseif(SelectFlag==1)
        CPre=FiltImg(CPre,6);
        CPre=FiltImg(CPre,4);
        CPre=imdilate(CPre,se);
    else
        C=FiltImg(C,6);
        C=FiltImg(C,4);
        C=imdilate(C,se);
    end
    CTotal=bitor(C,CPre);
    CTotal=bitor(CTotal,CPre2);
    CTotalTmp=CTotal;
    CTotal=bitand(CTotal,CTotalPre);
    CTotalPre=CTotalTmp;
    t=toc
    fprintf(text3,'%f\n',t);
    A=CellA{FrameIndex,1};
    FrameIndex=FrameIndex+1;
    if FrameIndex>TimeHeight
        FrameIndex=1;
    end
    FrameIndex1=FrameIndex1+1;
    if FrameIndex1>TimeHeight-1
        FrameIndex1=1;
    end
    SelectFlag=SelectFlag+1;
    if SelectFlag>2
        SelectFlag=0;
    end
%     title(num2str(f));
    C1(:,1:FrameWidth)=A;
    LeftX=fscanf(text1,'%f',1);
    TopY=fscanf(text1,'%f',1);
    RightX=fscanf(text1,'%f',1);
    BotY=fscanf(text1,'%f',1);
    GT=(RightX-LeftX+1)*(BotY-TopY+1);
    [nonzerorow nonzerocol]=find(CTotal>0);
    ForeTotal=length(nonzerorow);
    tmp=find(nonzerorow>=TopY&nonzerorow<=BotY&nonzerocol>=LeftX&nonzerocol<=RightX);
    TP=length(tmp);
    P=TP/ForeTotal;
    R=TP/GT;
    F1=2*TP/(GT+ForeTotal);
    fprintf(text2,'%f ',TP);
    fprintf(text2,'%f ',ForeTotal);
    fprintf(text2,'%f ',GT);
    fprintf(text2,'%f ',P);
    fprintf(text2,'%f ',R);
    fprintf(text2,'%f\n',F1);
    
    tmp=find(CTotal>0);
    C1(tmp)=255;
    Target(:,:,1)=C1;
    C1(tmp)=0;
    Target(:,:,2)=C1;
    Target(:,:,3)=C1;
    
%     figure(4);
%     imshow(CTotal);
%     figure(5);
%     imshow(Target);
%     title(num2str(f));

    %aviobj=addframe(aviobj,Target);
    %filename=sprintf('%s%d.jpg','resultimg\PTZ\',f);
    folder=dir('resultimg\SegCol');
    if isempty(folder)
        mkdir('resultimg\SegCol');
    end
    folder=dir('resultimg\SegBinary');
    if isempty(folder)
        mkdir('resultimg\SegBinary');
    end
    filename=sprintf('%s%d.bmp','resultimg\SegBinary\',f);
    imwrite(CTotal,filename);
    filename=sprintf('%s%d.bmp','resultimg\SegCol\',f);
    imwrite(Target,filename);
    f
end
fclose(text1);
fclose(text2);
fclose(text3);
%aviobj = close(aviobj);
a=0;
%       string='D:\Visual Tracking\Frequency moving target detection\resultimg\';
%       string=strcat(string,[num2str(f),'fftshift.bmp']);
%       imwrite(Tmp,string);

