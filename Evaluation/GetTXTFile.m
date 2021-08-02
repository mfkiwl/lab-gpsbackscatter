function [SaveDirName, FileTable] = GetTXTFile(SourceDirName)
%GETFILE 此处显示有关此函数的摘要
%   此处显示详细说明
W203Dir = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
% C12 TAG_1 10 序列
C12Dir_1 = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc1';
% C12 TAG_2 10 序列
C12Dir_2 = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc2';
% C12 TAG_2 PN 序列
C12Dir_2PN = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc2_PN';

if strcmp(SourceDirName, W203Dir)
    SaveDirName = './W203';
    FileTable = [
    'P01_150mV_100mV_Tag2_gnss_log_2021_06_12_17_13_53.txt';
    'P01_150mV_100mV_Tag2_gnss_log_2021_06_12_17_32_33.txt';
    'P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt';
    'P03_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_51.txt';
    'P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51.txt';
    'P05_150mV_100mV_Tag2_gnss_log_2021_06_12_18_46_52.txt';
    'P06_150mV_100mV_Tag2_gnss_log_2021_06_12_17_50_41.txt';
    % 'P06_150mv_100mV_Tag2_gnss_log_2021_06_12_18_59_13.txt';
    'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_46.txt';
    'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_19_10_10.txt';
    'P08_150mV_100mV_Tag2_gnss_log_2021_06_12_19_33_19.txt';
    'P08_150mV_100mV_Tag2_gnss_log_2021_06_12_19_58_04.txt';
    'P09_150mV_100mV_Tag2_gnss_log_2021_06_12_20_09_05.txt';
    'P09_150mV_100mV_Tag2_gnss_log_2021_06_12_20_54_12.txt';
    'P10_152mV_100mV_Tag2_gnss_log_2021_06_12_21_04_39.txt';
    'P12_150mV_100mV_Tag2_gnss_log_2021_06_12_21_15_30.txt';
    'P13_150mV_100mV_Tag2_gnss_log_2021_06_12_21_25_55.txt';
    'P14_150mV_100mV_Tag2_gnss_log_2021_06_12_21_36_49.txt';
    'P15_150mV_100mV_Tag2_gnss_log_2021_06_12_21_48_25.txt';
    'P16_150mV_100mV_Tag2_gnss_log_2021_06_12_22_10_28.txt';
    'P16_150mV_100mV_Tag2_gnss_log_2021_06_12_22_23_55.txt';
    'P17_150mV_100mV_Tag2_gnss_log_2021_06_12_22_34_52.txt';
    'P18_150mV_100mV_Tag2_gnss_log_2021_06_12_22_45_33.txt';
    'P19_150mV_100mV_Tag2_gnss_log_2021_06_12_22_56_08.txt';
    'P20_150mV_100mV_Tag2_gnss_log_2021_06_12_23_06_39.txt';
    'P21_150mV_100mV_Tag2_gnss_log_2021_06_13_20_17_14.txt';
    'P21_150mV_100mV_Tag2_gnss_log_2021_06_13_20_28_42.txt';
    'P22_150mV_100mV_Tag2_gnss_log_2021_06_13_20_39_21.txt';
    'P23_150mV_100mV_Tag2_gnss_log_2021_06_13_20_50_35.txt';
    'P24_150mV_100mV_Tag2_gnss_log_2021_06_13_21_01_53.txt';
    'P25_150mV_100mV_Tag2_gnss_log_2021_06_13_21_13_45.txt';
    'P26_150mV_100mV_Tag2_gnss_log_2021_06_13_21_25_00.txt';
    'P27_150mV_100mV_Tag2_gnss_log_2021_06_13_21_37_08.txt';
    'P28_150mV_100mV_Tag2_gnss_log_2021_06_13_21_49_03.txt';
    'P28_150mV_100mV_Tag2_gnss_log_2021_06_14_17_26_41.txt';
    'P29_150mV_100mV_Tag2_gnss_log_2021_06_13_22_09_41.txt';
    'P29_150mV_100mV_Tag2_gnss_log_2021_06_14_17_15_42.txt';
    'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_15_13_24.txt';
    'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_15_24_38.txt';
    'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_17_04_51.txt';
    % 'P31_150mV_100mV_Tag2_gnss_log_2021_06_14_15_42_07.txt'; 
    'P32_150mV_100mV_Tag2_gnss_log_2021_06_14_15_53_50.txt';
    'P33_150mV_100mV_Tag2_gnss_log_2021_06_14_16_31_26.txt';
    'P34_150mV_100mV_Tag2_gnss_log_2021_06_14_16_42_18.txt';
    'P35_150mV_100mV_Tag2_gnss_log_2021_06_14_16_54_20.txt';  
    ];

elseif strcmp(SourceDirName, C12Dir_1)
    SaveDirName = './C12_1';
    FileTable = [
        'P1R_150mV_xxmV_gnss_log_2021_07_28_09_52_03.txt';
        'P2R_150mV_xxmV_gnss_log_2021_07_28_09_59_16.txt';
        'P3R_150mV_xxmV_gnss_log_2021_07_28_12_58_17.txt';
        'P3R_150mV_xxmV_gnss_log_2021_07_28_13_13_51.txt';
        'P4R_150mV_xxmV_gnss_log_2021_07_28_10_22_52.txt';
        'P5R_150mV_xxmV_gnss_log_2021_07_28_10_33_26.txt';
        'P5R_150mV_xxmV_gnss_log_2021_07_28_10_46_16.txt';
        'P6R_150mV_xxmV_gnss_log_2021_07_28_11_08_33.txt';
        'P6R_150mV_xxmV_gnss_log_2021_07_28_11_18_51.txt';
        'P7R_150mV_xxmV_gnss_log_2021_07_28_11_29_44.txt';
        'P7R_150mV_xxmV_gnss_log_2021_07_28_11_40_54.txt';
        'P8R_150mV_xxmV_gnss_log_2021_07_28_11_52_38.txt';
        'P8R_150mV_xxmV_gnss_log_2021_07_28_12_03_18.txt';
        'P9R_150mV_xxmV_gnss_log_2021_07_28_12_21_53.txt';
        'P9R_150mV_xxmV_gnss_log_2021_07_28_12_32_51.txt';
    ];   
    
    elseif strcmp(SourceDirName, C12Dir_2)
    SaveDirName = './C12_2';
    FileTable = [
        'P1_150mV_xxmV_gnss_log_2021_07_28_15_09_37.txt';
        'P1_150mV_xxmV_gnss_log_2021_07_28_15_19_50.txt';
        'P2_150mV_xxmV_gnss_log_2021_07_28_15_33_04.txt';
        'P3_150mV_xxmV_gnss_log_2021_07_28_15_44_21.txt';
        'P4_150mV_xxmV_gnss_log_2021_07_28_15_55_50.txt';
        'P5_150mV_xxmV_gnss_log_2021_07_28_16_07_33.txt';
        'P6_150mV_xxmV_gnss_log_2021_07_28_16_18_05.txt';
        'P7_150mV_xxmV_gnss_log_2021_07_28_16_28_55.txt';
        'P8_150mV_xxmV_gnss_log_2021_07_28_16_39_57.txt';
        'P9_150mV_xxmV_gnss_log_2021_07_28_16_50_40.txt';
    ];    
    elseif strcmp(SourceDirName, C12Dir_2PN)
    SaveDirName = './C12_2PN';
    FileTable = [

    ]; 
else
    error('Directory has not been record !');
end

end

