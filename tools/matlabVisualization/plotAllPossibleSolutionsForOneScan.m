clc
clear



nameOfFolder = "out";

%%

set(groot,'defaultAxesTickLabelInterpreter','latex');  

magnitudeFFTW1 = readmatrix(nameOfFolder+"/magnitudeFFTW1.csv");
N=sqrt(size(magnitudeFFTW1,1));

phaseFFTW1 = readmatrix(nameOfFolder+"/phaseFFTW1.csv");
voxelDataUsed1 = readmatrix(nameOfFolder+"/voxelDataFFTW1.csv");

magnitudeFFTW2 = readmatrix(nameOfFolder+"/magnitudeFFTW2.csv");
phaseFFTW2 = readmatrix(nameOfFolder+"/phaseFFTW2.csv");
voxelDataUsed2 = readmatrix(nameOfFolder+"/voxelDataFFTW2.csv");
%fftwResult = complex(realFFTW,imaginaryFFTW);
%fftwResult = (reshape(fftwResult,N,N,N));
magnitude1 =zeros(N,N);
phase1 =zeros(N,N);
voxelData1 =zeros(N,N);
magnitude2 =zeros(N,N);
phase2 =zeros(N,N);
voxelData2 =zeros(N,N);
for j =1:N
    for k =1:N
        magnitude1(k,j) = magnitudeFFTW1(k*N-N+j);
        phase1(k,j) = phaseFFTW1(k*N-N+j);
        voxelData1(k,j) = voxelDataUsed1((k-1)*N+j);
        magnitude2(k,j) = magnitudeFFTW2(k*N-N+j);
        phase2(k,j) = phaseFFTW2(k*N-N+j);
        voxelData2(k,j) = voxelDataUsed2((k-1)*N+j);
    end
end


magnitude1=fftshift(magnitude1);
phase1=fftshift(phase1);
magnitude2=fftshift(magnitude2);
phase2=fftshift(phase2);

figure(1)
clf
% sigma = 70;
% size2 = size(magnitude1,1);
% [X, Y] = meshgrid(1:size2, 1:size2);
% G = 1*exp(-1/(sigma^2)*((Y-size2/2).^2 + (X-size2/2).^2));
% G = 1-G;
% imagesc(squeeze(G));

imagesc(squeeze(magnitude1));

%title('Magnitude first PCL Voxel:', 'Interpreter', 'latex')
axis image
box on 
pbaspect([1 1 1])

figure(2)
clf
imagesc((voxelData1))
axis image
%title('First PointCloud as Voxel:', 'Interpreter', 'latex')
% figure 2
pbaspect([1 1 1])

figure(3)
clf

imagesc(squeeze(magnitude2));


%imshow(magnitude2(:,:,size(magnitude2,1)/2));
%title('Magnitude second PCL Voxel:', 'Interpreter', 'latex')
axis image
pbaspect([1 1 1])

figure(4)
imagesc((voxelData2))
axis image
%title('Second PointCloud as Voxel:', 'Interpreter', 'latex')
pbaspect([1 1 1])
%%
resampledDataForSphere1 = readmatrix(nameOfFolder+"/resampledVoxel1.csv");
resampledDataForSphere2 = readmatrix(nameOfFolder+"/resampledVoxel2.csv");

resampledDataForSphereResult1 =zeros(N,N);
resampledDataForSphereResult2 =zeros(N,N);
for j = 1:N
    for i =1:N
            resampledDataForSphereResult1(j,i) = resampledDataForSphere1((i-1)*N+j);
            resampledDataForSphereResult2(j,i) = resampledDataForSphere2((i-1)*N+j);
    end
end

% currently not that interesting
if 1
    figure(5)
    clf

    sphere(1000);
    ch = get(gca,'children');
    set(ch,'facecolor','texturemap','cdata',resampledDataForSphereResult2,'edgecolor','none');
    axis equal
    %title('Second Resampled:', 'Interpreter', 'latex')
    pbaspect([1 1 1])
    figure(6)
    clf
    sphere(1000);
    ch = get(gca,'children');
    set(ch,'facecolor','texturemap','cdata',resampledDataForSphereResult1,'edgecolor','none');
    axis equal
    %title('first Resampled:', 'Interpreter', 'latex')
    pbaspect([1 1 1])
   
else
    figure(5)
    clf
    imagesc((resampledDataForSphereResult2));
    axis image
        %title('second Resampled:', 'Interpreter', 'latex')
    pbaspect([1 1 1])
    figure(6)
    clf
    imagesc((resampledDataForSphereResult1));
    axis image
    pbaspect([1 1 1])
        %title('first Resampled:', 'Interpreter', 'latex')
end

%%
figure(7)
correlationOfAngles = readmatrix(nameOfFolder+"/resultingCorrelation1D.csv");

anglesX = linspace(0,2*pi,size(correlationOfAngles,1));
[psor,lsor] = findpeaks(correlationOfAngles);
plot(anglesX,correlationOfAngles)
for i=1:size(psor,1)
    text(anglesX(lsor(i))+0.15,psor(i)+1,[num2str(round(anglesX(lsor(i)),2)) ], 'Interpreter', 'latex')
end
grid on
xlim([-0.2 2*pi+0.2])
%ylim([min(correlationOfAngles)-10 max(correlationOfAngles)+10])
%title('Correlation of different Angles', 'Interpreter', 'latex')
box on 
xlabel("angle in rad", 'Interpreter', 'latex');
ylabel("Correlation Value", 'Interpreter', 'latex');
pbaspect([1 1 1])
%%




%%

dataInformation = readmatrix(nameOfFolder+"/dataForReadIn.csv");

numberOfSolutions = dataInformation(1);
bestSolution = dataInformation(2);





% [M,I] = max(correlationMatrixShift2D, [], "all", "linear");
% [dim1, dim2, dim3] = ind2sub(size(correlationMatrixShift2D),I)

%correlationMatrixShift2D  = squeeze(correlationMatrixShift3D(65,:,:));
figure(8)
clf

for i = 1:numberOfSolutions
        subplot(2,2,i)
    correlationMatrixShift1D = readmatrix(nameOfFolder+"/resultingCorrelationShift"+ string(num2str(i-1)) +".csv");
    resultSize = nthroot(length(correlationMatrixShift1D),2);
    %A = reshape(results,resultSize,resultSize);
    correlationMatrixShift2D = reshape(correlationMatrixShift1D,resultSize,resultSize);
%     imagesc(correlationMatrixShift2D);
    [Xplot,Yplot]=meshgrid(1:resultSize,1:resultSize);
    surf(Xplot,Yplot,(correlationMatrixShift2D),'edgecolor', 'none');
    %xlabel("x-axis", 'Interpreter', 'latex');
    %ylabel("y-axis", 'Interpreter', 'latex');
    
    %title('Correlation Surface of  Best Performing Match', 'Interpreter', 'latex')
    box on 
    title(num2str(i), 'Interpreter', 'latex')
    pbaspect([1 1 1])
end
%set(gcf,'units','points','position',[0,0,500,500])





%%



f = figure(9);
clf




for i = 1:numberOfSolutions
    resultVoxel1TMP= readmatrix(nameOfFolder+"/resultVoxel1"+num2str(i-1)+".csv");
    resultVoxel2TMP= readmatrix(nameOfFolder+"/resultVoxel2"+num2str(i-1)+".csv");

    %fftwResult = complex(realFFTW,imaginaryFFTW);
    %fftwResult = (reshape(fftwResult,N,N,N));
    voxelResult1 =zeros(N,N);
    voxelResult2 =zeros(N,N);
    for j =1:N
        for k =1:N
            voxelResult1(k,j) = resultVoxel1TMP(k*N-N+j);
            voxelResult2(k,j) = resultVoxel2TMP(k*N-N+j);
    
        end
    end



    subplot(2,2,i)
    C = imfuse(voxelResult1,voxelResult2,'blend','Scaling','joint');
    imagesc(squeeze(C));

    title(num2str(i), 'Interpreter', 'latex')

    box on


end



