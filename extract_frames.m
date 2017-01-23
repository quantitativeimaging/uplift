% Read frames at user-specified times in MP4 input and write as png

filenameMP4 = [pwd,'\input\Nov 22 Exp 7 H180 D45.MP4']; % Name of file to process

tInit = 12.75; % Seconds. Time value of initial frame
tPeak = 47.4; % Seconds. Time value of peak force frame

v = VideoReader(filenameMP4);
v.CurrentTime = tInit;
imDat = readFrame(v);
imwrite(imDat, [pwd, '\output_initial_and_peak_frames\frame_init.png'] );

v.CurrentTime = tPeak;
imDat = readFrame(v);
imwrite(imDat, [pwd, '\output_initial_and_peak_frames\frame_peak.png'] );

