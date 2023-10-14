datapath = '/Users/qihuilyu/Desktop/Data/HDRCT/T2';
ctpath = fullfile(datapath,'CT');
RTdosepath = fullfile(datapath,'RTDOSE');
RTplanpath = fullfile(datapath,'RTPLAN');
RTstpath = fullfile(datapath,'RTst');


[V,spatial,dim] = dicomreadVolume(ctpath);
V = squeeze(V);
figure;imshow3D(flip(V,1),[800,1500])

RTstfile = dir(fullfile(RTstpath,'**/*.dcm'));
file = fullfile(RTstpath, RTstfile.name);
infoRTst = dicominfo(file,'UseVRHeuristic',false);
rtContours = dicomContours(infoRTst);
plotContour(rtContours)

referenceInfo = imref3d(size(V),xlim,ylim,zlim);
contourIndex = 21;  
rtMask = createMask(rtContours, contourIndex, spatial);
figure;imshow3D(rtMask,[]);

figure;imshow3D([rtMask, double(V)/double(max(V(:)))],[]);



CTfiles = dir(fullfile(ctpath,'**/*.dcm'));
[sortedFile,index] = sort_nat({CTfiles.name});

for ii = 1:numel(ctpath)
    file = fullfile(ctpath, CTfiles(ii).name);
    CT(:,:,ii) = dicomread(file);
    infoCT{ii} = dicominfo(file);

    DICOMImageStruct(ii) = infoCT{ii};
    DICOMimagefilenames{ii} = CTfiles(ii).name;
end

RTdosefiles = dir(RTdosepath);
file = fullfile(RTdosepath, RTdosefiles(3).name);
RTdose = dicomread(file);
infoRTdose = dicominfo(file);

RTplanfiles = dir(RTplanpath);
file = fullfile(RTplanpath, RTplanfiles(3).name);
RTplan = dicomread(file);
infoRTplan = dicominfo(file);

RTstfiles = dir(RTstpath);
file = fullfile(RTstpath, RTstfiles(3).name);
RTst = dicomread(file);
infoRTst = dicominfo(file);


         RTstructfiles = FileName;
            contour0 = dicomContours(infoRTst);
            referenceInfo = imref3d([512,512,40],xlim,ylim,zlim);
            referenceInfo = imref3d;
            contourIndex = 1;  
            rtMask = createMask(contour0, contourIndex, referenceInfo);
            






structures = LoadDICOMStructures(DicomPath,RTstructfiles,image);

    count = 1;
    count2 = 1;
    for ii = 1:length(baseFileNames)
        FileName = baseFileNames(index(ii)).name;
        fullFileName = fullfile(DicomPath, FileName);
        info = dicominfo(fullFileName);
        if(strcmp(info.Modality,'CT'))
            DICOMImageStruct(ii) = info;
            DICOMimagefilenames{count} = FileName;
            count = count+1;
        end
        if(strcmp(info.Modality,'RTSTRUCT'))
            RTstructfiles = FileName;
            contour0 = dicomContours(info);
            referenceInfo = imref3d([128,128,50],xlim,ylim,zlim);
            referenceInfo = imref3d;
            contourIndex = 1;  
            rtMask = createMask(contour0, contourIndex, referenceInfo);
            
            if(count2>1)
                error('Multiple RT structure files!')
            end
            count2 = count2 + 1;
        end
    end
    
    image = LoadDICOMImages(DicomPath,...
        DICOMimagefilenames);
    
