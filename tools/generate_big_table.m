
   for s = 1:length(coh) % das hier über coh machen
    % subject ID
    clear patient_id
    patient_id = coh{s,1}.id; % und hier dann die ID 
    
    % create empty cell
    cell_info = cell(1,129); % one person, 129 values/info aspects
    
for row = 1
   for col = 1
      cell_info{row,col} = patient_id;% VPCode
   end 
   for col = 2
       cell_info{row,col} = coh_delta{s};% mean coherence delta
   end
   for col = 3
       cell_info{row,col} = coh_beta{s};% mean coherence beta
   end
   for col = 4
       cell_info{row,col} = threshold_delta_01{s};% threshold delta 0.1
   end
   for col = 5
       cell_info{row,col} = threshold_delta_02{s};% threshold delta 0.2
   end
   for col = 6
       cell_info{row,col} = threshold_delta_03{s};%
   end
   for col = 7
       cell_info{row,col} = threshold_delta_04{s};% 
   end
   for col = 8
       cell_info{row,col} = threshold_delta_05{s};% 
   end 
   for col = 9
       cell_info{row,col} = threshold_delta_06{s};%
   end
   for col = 10
       cell_info{row,col} = threshold_delta_07{s};% 
   end
   for col = 11
       cell_info{row,col} = threshold_delta_08{s};% 
   end
   for col = 12
       cell_info{row,col} = threshold_delta_09{s};% 
   end
   for col = 13
       cell_info{row,col} = threshold_beta_01{s};% 
   end
   for col = 14
       cell_info{row,col} = threshold_beta_02{s};% 
   end
   for col = 15
       cell_info{row,col} = threshold_beta_03{s};% 
   end
   for col = 16
       cell_info{row,col} = threshold_beta_04{s};% 
   end
   for col = 17
       cell_info{row,col} = threshold_beta_05{s};% 
   end
   for col = 18
       cell_info{row,col} = threshold_beta_06{s};% 
   end
   for col = 19
       cell_info{row,col} = threshold_beta_07{s};% 
   end
   for col = 20
       cell_info{row,col} = threshold_beta_08{s};% 
   end
   for col = 21
       cell_info{row,col} = threshold_beta_09{s};% 
   end
   for col = 22
       cell_info{row,col} = gcc_delta_01{s};%
   end
   for col = 23
       cell_info{row,col} = gcc_delta_02{s};% 
   end
   for col = 24
       cell_info{row,col} = gcc_delta_03{s};% 
   end
   for col = 25
       cell_info{row,col} = gcc_delta_04{s};% 
   end
   for col = 26
       cell_info{row,col} = gcc_delta_05{s};% 
   end
   for col = 27
       cell_info{row,col} = gcc_delta_06{s};% 
   end
   for col = 28
       cell_info{row,col} = gcc_delta_07{s};% 
   end
   for col = 29
       cell_info{row,col} = gcc_delta_08{s};% 
   end
   for col = 30
       cell_info{row,col} = gcc_delta_09{s};% 
   end
   for col = 31
       cell_info{row,col} = gcc_beta_01{s};% 
   end
   for col = 32
       cell_info{row,col} = gcc_beta_02{s};%
   end
   for col = 33
       cell_info{row,col} = gcc_beta_03{s};% 
   end
   for col = 34
       cell_info{row,col} = gcc_beta_04{s};% 
   end
   for col = 35
       cell_info{row,col} = gcc_beta_05{s};% 
   end
   for col = 36
       cell_info{row,col} = gcc_beta_06{s};%
   end
   for col = 37
       cell_info{row,col} = gcc_beta_07{s};% 
   end
   for col = 38
       cell_info{row,col} = gcc_beta_08{s};% 
   end
   for col = 39
       cell_info{row,col} = gcc_beta_09{s};% 
   end
   for col = 40
       cell_info{row,col} = cpl_delta_01{s};% 
   end
   for col = 41
       cell_info{row,col} = cpl_delta_02{s};% 
   end
   for col = 42
       cell_info{row,col} = cpl_delta_03{s};% 
   end
   for col = 43
       cell_info{row,col} = cpl_delta_04{s};% 
   end
   for col = 44
       cell_info{row,col} = cpl_delta_05{s};% 
   end
   for col = 45
       cell_info{row,col} = cpl_delta_06{s};% 
   end
   for col = 46
       cell_info{row,col} = cpl_delta_07{s};% 
   end
   for col = 47
       cell_info{row,col} = cpl_delta_08{s};% 
   end
   for col = 48
       cell_info{row,col} = cpl_delta_09{s};% 
   end
   for col = 49
       cell_info{row,col} = cpl_beta_01{s};% 
   end
   for col = 50
       cell_info{row,col} = cpl_beta_02{s};% 
   end
   for col = 51
       cell_info{row,col} = cpl_beta_03{s};% 
   end
   for col = 52
       cell_info{row,col} = cpl_beta_04{s};% 
   end
   for col = 53
       cell_info{row,col} = cpl_beta_05{s};% 
   end
   for col = 54
       cell_info{row,col} = cpl_beta_06{s};% 
   end
   for col = 55
       cell_info{row,col} = cpl_beta_07{s};% 
   end
   for col = 56
       cell_info{row,col} = cpl_beta_08{s};% 
   end
   for col = 57
       cell_info{row,col} = cpl_beta_09{s};% 
   end
   for col = 58
       cell_info{row,col} = gcc_rand_delta_01{s};% 
   end
   for col = 59
       cell_info{row,col} = gcc_rand_delta_02{s};% 
   end
   for col = 60
       cell_info{row,col} = gcc_rand_delta_03{s};% 
   end
   for col = 61
       cell_info{row,col} = gcc_rand_delta_04{s};% 
   end
   for col = 62
       cell_info{row,col} = gcc_rand_delta_05{s};% 
   end
   for col = 63
       cell_info{row,col} = gcc_rand_delta_06{s};% 
   end
   for col = 64
       cell_info{row,col} = gcc_rand_delta_07{s};% 
   end
   for col = 65
       cell_info{row,col} = gcc_rand_delta_08{s};% 
   end
   for col = 66
       cell_info{row,col} = gcc_rand_delta_09{s};% 
   end
   for col = 67
       cell_info{row,col} = gcc_rand_beta_01{s};% 
   end
   for col = 68
       cell_info{row,col} = gcc_rand_beta_02{s};% 
   end
   for col = 69
       cell_info{row,col} = gcc_rand_beta_03{s};% 
   end
   for col = 70
       cell_info{row,col} = gcc_rand_beta_04{s};% 
   end
   for col = 71
       cell_info{row,col} = gcc_rand_beta_05{s};% 
   end
   for col = 72
       cell_info{row,col} = gcc_rand_beta_06{s};% 
   end
   for col = 73
       cell_info{row,col} = gcc_rand_beta_07{s};% 
   end
   for col = 74
       cell_info{row,col} = gcc_rand_beta_08{s};% 
   end
   for col = 75
       cell_info{row,col} = gcc_rand_beta_09{s};% 
   end
   for col = 76
       cell_info{row,col} = cpl_rand_delta_01{s};% 
   end
   for col = 77
       cell_info{row,col} = cpl_rand_delta_02{s};% 
   end
   for col = 78
       cell_info{row,col} = cpl_rand_delta_03{s};% 
   end
   for col = 79
       cell_info{row,col} = cpl_rand_delta_04{s};% 
   end
   for col = 80
       cell_info{row,col} = cpl_rand_delta_05{s};% 
   end
   for col = 81
       cell_info{row,col} = cpl_rand_delta_06{s};% 
   end
   for col = 82
       cell_info{row,col} = cpl_rand_delta_07{s};% 
   end
   for col = 83
       cell_info{row,col} = cpl_rand_delta_08{s};% 
   end
   for col = 84
       cell_info{row,col} = cpl_rand_delta_09{s};% 
   end
   for col = 85
       cell_info{row,col} = cpl_rand_beta_01{s};% 
   end
   for col = 86
       cell_info{row,col} = cpl_rand_beta_02{s};% 
   end
   for col = 87
       cell_info{row,col} = cpl_rand_beta_03{s};% 
   end
   for col = 88
       cell_info{row,col} = cpl_rand_beta_04{s};% 
   end
   for col = 89
       cell_info{row,col} = cpl_rand_beta_05{s};% 
   end
   for col = 90
       cell_info{row,col} = cpl_rand_beta_06{s};% 
   end
   for col = 91
       cell_info{row,col} = cpl_rand_beta_07{s};% 
   end
   for col = 92
       cell_info{row,col} = cpl_rand_beta_08{s};% 
   end
   for col = 93
       cell_info{row,col} = cpl_rand_beta_09{s};% 
   end
   for col = 94
       cell_info{row,col} = smallworldness_delta_01{s};% 
   end
   for col = 95
       cell_info{row,col} = smallworldness_delta_02{s};% 
   end
   for col = 96
       cell_info{row,col} = smallworldness_delta_03{s};% 
   end
   for col = 97
       cell_info{row,col} = smallworldness_delta_04{s};% 
   end
   for col = 98
       cell_info{row,col} = smallworldness_delta_05{s};% 
   end
   for col = 99
       cell_info{row,col} = smallworldness_delta_06{s};% 
   end
   for col = 100
       cell_info{row,col} = smallworldness_delta_07{s};% 
   end
   for col = 101
       cell_info{row,col} = smallworldness_delta_08{s};% 
   end
   for col = 102
       cell_info{row,col} = smallworldness_delta_09{s};% 
   end
   for col = 103
       cell_info{row,col} = smallworldness_beta_01{s};% 
   end
   for col = 104
       cell_info{row,col} = smallworldness_beta_02{s};% 
   end
   for col = 105
       cell_info{row,col} = smallworldness_beta_03{s};% 
   end
   for col = 106
       cell_info{row,col} = smallworldness_beta_04{s};% 
   end
   for col = 107
       cell_info{row,col} = smallworldness_beta_05{s};% 
   end
   for col = 108
       cell_info{row,col} = smallworldness_beta_06{s};% 
   end
   for col = 109
       cell_info{row,col} = smallworldness_beta_07{s};% 
   end
   for col = 110
       cell_info{row,col} = smallworldness_beta_08{s};% 
   end
   for col = 111
       cell_info{row,col} = smallworldness_beta_09{s};% 
   end
   for col = 112
       cell_info{row,col} = coh_delta_01{s};
   end
   for col = 113
       cell_info{row,col} = coh_delta_02{s};
   end
   for col = 114
       cell_info{row,col} = coh_delta_03{s};
   end
   for col = 115
       cell_info{row,col} = coh_delta_04{s};
   end
   for col = 116
       cell_info{row,col} = coh_delta_05{s};
   end
   for col = 117
       cell_info{row,col} = coh_delta_06{s};
   end
   for col = 118
       cell_info{row,col} = coh_delta_07{s};
   end
   for col = 119
       cell_info{row,col} = coh_delta_08{s};
   end
   for col = 120
       cell_info{row,col} = coh_delta_09{s};
   end
   for col = 121
       cell_info{row,col} = coh_beta_01{s};
   end
   for col = 122
       cell_info{row,col} = coh_beta_02{s};
   end
   for col = 123
       cell_info{row,col} = coh_beta_03{s};
   end
   for col = 124
       cell_info{row,col} = coh_beta_04{s};
   end
   for col = 125
       cell_info{row,col} = coh_beta_05{s};
   end
   for col = 126
       cell_info{row,col} = coh_beta_06{s};
   end
   for col = 127
       cell_info{row,col} = coh_beta_07{s};
   end
   for col = 128
       cell_info{row,col} = coh_beta_08{s};
   end
   for col = 129
       cell_info{row,col} = coh_beta_09{s};
   end
end 

% create the table names
VarNames = ["participant_id" "coh_delta" "coh_beta" "thr_delta_01" "thr_delta_02" "thr_delta_03" "thr_delta_04" "thr_delta_05" "thr_delta_06" "thr_delta_07" "thr_delta_08" "thr_delta_09" "thr_beta_01" "thr_beta_02" "thr_beta_03" "thr_beta_04" "thr_beta_05" "thr_beta_06" "thr_beta_07" "thr_beta_08" "thr_beta_09" "gcc_delta_01" "gcc_delta_02" "gcc_delta_03" "gcc_delta_04" "gcc_delta_05" "gcc_delta_06" "gcc_delta_07" "gcc_delta_08" "gcc_delta_09" "gcc_beta_01" "gcc_beta_02" "gcc_beta_03" "gcc_beta_04" "gcc_beta_05" "gcc_beta_06" "gcc_beta_07" "gcc_beta_08" "gcc_beta_09" "cpl_delta_01" "cpl_delta_02" "cpl_delta_03" "cpl_delta_04" "cpl_delta_05" "cpl_delta_06" "cpl_delta_07" "cpl_delta_08" "cpl_delta_09" "cpl_beta_01" "cpl_beta_02" "cpl_beta_03" "cpl_beta_04" "cpl_beta_05" "cpl_beta_06" "cpl_beta_07" "cpl_beta_08" "cpl_beta_09" "gcc_rand_delta_01" "gcc_rand_delta_02" "gcc_rand_delta_03" "gcc_rand_delta_04" "gcc_rand_delta_05" "gcc_rand_delta_06" "gcc_rand_delta_07" "gcc_rand_delta_08" "gcc_rand_delta_09" "gcc_rand_beta_01" "gcc_rand_beta_02" "gcc_rand_beta_03" "gcc_rand_beta_04" "gcc_rand_beta_05" "gcc_rand_beta_06" "gcc_rand_beta_07" "gcc_rand_beta_08" "gcc_rand_beta_09" "cpl_rand_delta_01" "cpl_rand_delta_02" "cpl_rand_delta_03" "cpl_rand_delta_04" "cpl_rand_delta_05" "cpl_rand_delta_06" "cpl_rand_delta_07" "cpl_rand_delta_08" "cpl_rand_delta_09" "cpl_rand_beta_01" "cpl_rand_beta_02" "cpl_rand_beta_03" "cpl_rand_beta_04" "cpl_rand_beta_05" "cpl_rand_beta_06" "cpl_rand_beta_07" "cpl_rand_beta_08" "cpl_rand_beta_09" "sw_delta_01" "sw_delta_02" "sw_delta_03" "sw_delta_04" "sw_delta_05" "sw_delta_06" "sw_delta_07" "sw_delta_08" "sw_delta_09" "sw_beta_01" "sw_beta_02" "sw_beta_03" "sw_beta_04" "sw_beta_05" "sw_beta_06" "sw_beta_07" "sw_beta_08" "sw_beta_09" "coh_delta_01" "coh_delta_02" "coh_delta_03" "coh_delta_04" "coh_delta_05" "coh_delta_06" "coh_delta_07" "coh_delta_08" "coh_delta_09" "coh_beta_01" "coh_beta_02" "coh_beta_03" "coh_beta_04" "coh_beta_05" "coh_beta_06" "coh_beta_07" "coh_beta_08" "coh_beta_09"];% muss noch vervollständigt werden

% gather data in a temporary table
currTable = table(cell_info(:,1),cell_info(:,2),cell_info(:,3),cell_info(:,4),cell_info(:,5),cell_info(:,6),cell_info(:,7),cell_info(:,8),cell_info(:,9),cell_info(:,10),cell_info(:,11),cell_info(:,12),cell_info(:,13),cell_info(:,14),cell_info(:,15),cell_info(:,16),cell_info(:,17),cell_info(:,18),cell_info(:,19),cell_info(:,20),cell_info(:,21),cell_info(:,22),cell_info(:,23),cell_info(:,24),cell_info(:,25),cell_info(:,26),cell_info(:,27),cell_info(:,28),cell_info(:,29),cell_info(:,30),cell_info(:,31),cell_info(:,32),cell_info(:,33),cell_info(:,34),cell_info(:,35),cell_info(:,36),cell_info(:,37),cell_info(:,38),cell_info(:,39),cell_info(:,40),cell_info(:,41),cell_info(:,42),cell_info(:,43),cell_info(:,44),cell_info(:,45),cell_info(:,46),cell_info(:,47),cell_info(:,48),cell_info(:,49),cell_info(:,50),cell_info(:,51),cell_info(:,52),cell_info(:,53),cell_info(:,54),cell_info(:,55),cell_info(:,56),cell_info(:,57),cell_info(:,58),cell_info(:,59),cell_info(:,60),cell_info(:,61),cell_info(:,62),cell_info(:,63),cell_info(:,64),cell_info(:,65),cell_info(:,66),cell_info(:,67),cell_info(:,68),cell_info(:,69),cell_info(:,70),cell_info(:,71),cell_info(:,72),cell_info(:,73),cell_info(:,74),cell_info(:,75),cell_info(:,76),cell_info(:,77),cell_info(:,78),cell_info(:,79),cell_info(:,80),cell_info(:,81),cell_info(:,82),cell_info(:,83),cell_info(:,84),cell_info(:,85),cell_info(:,86),cell_info(:,87),cell_info(:,88),cell_info(:,89),cell_info(:,90),cell_info(:,91),cell_info(:,92),cell_info(:,93),cell_info(:,94),cell_info(:,95),cell_info(:,96),cell_info(:,97),cell_info(:,98),cell_info(:,99),cell_info(:,100),cell_info(:,101),cell_info(:,102),cell_info(:,103),cell_info(:,104),cell_info(:,105),cell_info(:,106),cell_info(:,107),cell_info(:,108),cell_info(:,109),cell_info(:,110),cell_info(:,111),cell_info(:,112),cell_info(:,113),cell_info(:,114),cell_info(:,115),cell_info(:,116),cell_info(:,117),cell_info(:,118),cell_info(:,119),cell_info(:,120),cell_info(:,121),cell_info(:,122),cell_info(:,123),cell_info(:,124),cell_info(:,125),cell_info(:,126),cell_info(:,127),cell_info(:,128),cell_info(:,129),'VariableNames',VarNames);
    
% Vertically concatenate the new table to the big table
bigTable = vertcat(bigTable, currTable);

   end