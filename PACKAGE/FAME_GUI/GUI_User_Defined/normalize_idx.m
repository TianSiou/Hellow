function [S]=normalize_idx(S1)
%% 規範化三圍做標
        %process string from answer( sphere center)
        S1  = strrep(S1, ',',' ');
        S1  = strrep(S1,']','');
        S1  = strrep(S1,'[','');
        [ t1, t2 ] = strtok(S1);
        [ t2, t3 ] = strtok(t2);
              t3   = strrep(t3,' ','');
        S  = [t1,' ',t2,' ',t3];
end
