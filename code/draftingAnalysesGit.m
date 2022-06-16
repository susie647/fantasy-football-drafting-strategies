%% drafting analyses for sleeperbot data


clear; close all;
doPrint = false;

% Analysis List
aList = {...
    'valueCliffs'                      ; ...
    % 'visualize'                   ; ...
    %   'playersInset'                          ; ...
    %  'roleOrder'                           ; ...
    %    'positionChoice'                    ; ...
    %    'cumulativePicks'                  ; ...
    % 'composition'                      ; ...
    %   'copycat3'                          ; ...
    %   'fandom2'                            ; ...
    %    'positionPerformance'                ; ...
    %'handcuff'                         ; ...
    %'descriptiveStats'                   ; ...
    %  'humansVsBots'                       ; ...
    %    'copycatConsequences'                ; ...
    };

%% constants
figDir = 'figures/';
dataDir = '../data/';

%% data
T = readtable([dataDir 'cleaned_data.csv'], 'headerlines', 1);

[n, ~] = size(T);

d.league = T{:, 1};
d.nLeagues = length(unique(d.league));
d.pick = T{:, 2};
d.team = T{:, 3};
d.leagueTeams = T{:, 5};
d.leagueHumans = T{:, 6};
d.player = T{:, 7};
d.playerNames = T{:, 9};
d.nPlayers = length(unique(d.player));

d.auto = zeros(n, 1);
d.auto(strmatch('True', T{:, 4})) = 1;

d.role = nan(n, 1);
d.roleBrief = nan(n, 1);
d.roleNames = {'QB', 'RB', 'WR', 'TE', 'K', 'DS'};
d.roleNamesBrief =  {'QB', 'RW', 'TE', 'K', 'DS'};
d.nRoles = numel(d.roleNames);
d.nRolesBrief = numel(d.roleNamesBrief);
for idx = 1:numel(d.roleNames)
    d.role(strmatch(d.roleNames{idx}, T{:, 8})) = idx;
end
for idx = 1:numel(d.roleNamesBrief)
    d.roleBrief(strmatch(d.roleNamesBrief{idx}, T{:, 8})) = idx;
end

d.wins = T{:, 11};
d.games = T{:, 12};

nflTeamRaw = T{:, 10};
d.nflTeamNames = unique(nflTeamRaw);
d.nflTeamNames = d.nflTeamNames(2:end);
d.nflTeam = nan(n, 1);
for i = 1:n
    match = find(strcmp(nflTeamRaw{i}, d.nflTeamNames));
    if ~isempty(match)
        d.nflTeam(i) = match;
    end
end

%% colors
load pantoneColors;

roleColorsFull = { ...
    pantone.Marsala        ; ... % QB
    pantone.ClassicBlue    ; ... % RB
    pantone.Treetop        ; ... % WR
    pantone.DuskBlue       ; ... % TE
    pantone.Titanium       ; ....% K
    pantone.Sandstone      ;     % DEF
    };


roleColorsBrief = { ...
    pantone.Marsala        ; ... % QB
    pantone.Comfrey     ; ... % RW
    pantone.DuskBlue       ; ... % TE
    pantone.Titanium       ; ....% K
    pantone.Sandstone      ;     % DEF
    };

winsColor = pantone.ClassicBlue;
roleMarkers = {'s', '^', 'v', 's', '^', 'v'};

%% loop over analyses
for aListi = 1:numel(aList)
    
    figAppend = '';
    
    switch aList{aListi}
        
        case 'composition'
            
            sizeTeamList = 15;
            collapseRBWR = true;
            proportionInclude = 0.95; % how many compositions

            fontSize = 16;
            w = 0.8; h = 0.8;
            maxConfigs = 1e4;
            
            if collapseRBWR
                nRoles = d.nRolesBrief;
                roleNames = d.roleNamesBrief;
                bins = [0.5 1.5 3.5 4.5 5.5 6.5];
                roleColors = roleColorsBrief;
            else
                nRoles = d.nRoles;
                roleNames = d.roleNames;
                bins = 0.5:6.5;
                roleColors = roleColorsFull;
            end
            
            leagueList = unique(d.league);
            
            for sizeIdx = 1:length(sizeTeamList)
                sizeTeam = sizeTeamList(sizeIdx);
                
                nCompositions = 0;
                m = nan(maxConfigs, nRoles);
                mCount = zeros(maxConfigs, 1);
                wins = cell(maxConfigs, 1);
                for leagueIdx = 1:length(leagueList)
                    match = find(d.league == leagueList(leagueIdx));
                    teamList = unique(d.team(match));
                    teamList = teamList(find(teamList >= 0));
                    for teamIdx = 1:length(teamList)
                        match = find(d.league == leagueList(leagueIdx) & d.team == teamList(teamIdx));
                        if ~isempty(match)
                            go = true;
                            if length(match) ~= sizeTeam
                                go = false;
                            end
                            if go
                                count = histcounts(d.role(match), bins);
                                if sum(count) == sizeTeam
                                    if leagueIdx == 1
                                        nCompositions = nCompositions + 1;
                                        m(1, :) = count; mCount(1) = 1;
                                        wins{1} = max(d.wins(match))/max(d.games(match)); % max needed because autopicked entries don't record wins (look row 187)
                                    else
                                        [isThere, matchIdx] = ismember(count, m, 'rows');
                                        if isThere
                                            mCount(matchIdx) = mCount(matchIdx) + 1;
                                            wins{matchIdx} = [wins{matchIdx} max(d.wins(match))/max(d.games(match))];
                                            
                                            
                                        else
                                            nCompositions = nCompositions + 1;
                                            m(nCompositions, :) = count; mCount(nCompositions) = 1;
                                            wins{nCompositions} = max(d.wins(match))/max(d.games(match));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                [~, sortIdx] = sort(mCount, 'descend');
                m = m(sortIdx(1:nCompositions), :);
                mCount = mCount(sortIdx(1:nCompositions));
                wins = wins(sortIdx(1:nCompositions));
                
                tmpProp = cumsum(mCount./sum(mCount));
                nCompositions = sum(tmpProp <= proportionInclude);
                
                % figure window
                F = figure; clf; hold on;
                setFigure(F, [0.2 0.1 0.55 0.8], '');
                
                % composition axes
                set(gca, ...
                    'units'         ,  'normalized'           , ...
                    'position'      ,  [0.1 0.075 0.8 0.65]   ,...
                    'xlim'          ,  [1/2 nCompositions+1/2]     , ...
                    'xtick'         ,  1:nCompositions            , ...
                    'xticklabel'    ,  [], ...
                    'ylim'          ,  [0 sizeTeam+3/2]       , ...
                    'ydir'          ,  'reverse'              , ...
                    'ycolor'        , 'none'          , ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.005 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'fontsize'      ,  9                      );
                xlabel('Team Composition', 'fontsize', fontSize+2);
                
                % plot
                for compIdx = 1:nCompositions
                    count = 1;
                    for roleIdx = 1:nRoles
                        if m(compIdx, roleIdx) > 0
                            for idx = 1:m(compIdx, roleIdx)
                                count = count + 1;
                                rectangle('position', [compIdx-w/2 count-h/2 w h], ...
                                    'facecolor' , [roleColors{roleIdx} 0.65] , ...
                                    'edgecolor' , 'w'                                             , ...
                                    'curvature' , [0.5 0.5]);
                                if idx == 1
                                    text(compIdx, count, lower(roleNames{roleIdx}), ...
                                        'fontsize', fontSize-4, ...
                                        'hor', 'cen', ...
                                        'vert', 'mid');
                                end
                            end
                        end
                    end
                    
                end
                
                % proportion axes
                pAX = axes; hold on;
                set(pAX, ...
                    'units'         ,  'normalized'           , ...
                    'position'      ,  [0.1 0.675 0.8 0.3]   ,...
                    'xlim'          ,  [1/2 nCompositions+1/2]     , ...
                    'xtick'         ,  1:nCompositions            , ...
                    'xticklabel'    , [], ...
                    'ylim'          ,  [0 0.3]       , ...
                    'ytick', 0:0.1:0.3, ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.005 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'fontsize'      ,  fontSize                      );
                
                H = bar(1:nCompositions,  mCount(1:nCompositions)/sum(mCount),  0.6, 'k');
                set(H, ...
                    'facecolor', pantone.Custard, ...
                    'edgecolor', 'none');
                ylabel('Proportion of Teams');
                
                wAX = axes; hold on;
                set(wAX, ...
                    'units'         ,  'normalized'           , ...
                    'position'      ,  [0.1 0.675 0.8 0.3]   ,...
                    'xlim'          ,  [1/2 nCompositions+1/2]     , ...
                    'xtick'         ,  []            , ...
                    'xticklabel'    , [], ...
                    'ylim'          ,  [0 1]       , ...
                    'ytick', 0:0.5:1, ...
                    'yaxisloc', 'right', ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.005 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'color'    , 'none', ...
                    'fontsize'      ,  fontSize                      );
                
                plot([0 nCompositions+1], [0.5 0.5], '--', ...
                    'color', pantone.Titanium, ...
                    'linewidth', 1);
                
                meanWins = cellfun(@nanmean, wins);
                seWins = cellfun(@nanstd, wins)./sqrt(cellfun(@length, wins));
                errorbar(1:nCompositions, meanWins(1:nCompositions), seWins(1:nCompositions), ...
                    'marker', 'o', ...
                    'markersize', 8, ...
                    'markerfacecolor', winsColor, ...
                    'markeredgecolor', 'w', ...
                    'color', winsColor);
                ylabel('Winning Proportion', 'rot', -90, 'vert', 'bot');
                
                
                fprintf('----\nTeams with %d players.\n', sizeTeam);
                if collapseRBWR
                    fprintf('RBs and WRs collapsed.\n');
                    figAppend = 'collapse';
                end
                fprintf('%d configurations needed to account for %d%%\n', nCompositions, proportionInclude*100);
                
                if doPrint
                    print(sprintf('figures/%s_%d%s.png', aList{aListi}, sizeTeam, figAppend), '-dpng', '-r600');
                    print(sprintf('%s/%s_%d%s.eps', figDir, aList{aListi}, sizeTeam, figAppend), '-depsc');
                end
                
            end
            
        case   'visualize'
            
            % user constants
            for whichLeague = 35:35
                
                fontSize = 12;
                width = 0.8; height = 0.9; gap = 0.25;
                
                uLeagues = unique(d.league);
                match = find(d.league == uLeagues(whichLeague));
                
                % derived constants
                nTeams = d.leagueTeams(match(1));
                picksPerTeam = length(match)/nTeams;
                teamList = unique(d.team(match));
                %             teamList = teamList(find(teamList >= 0));
                sortIdx = d.team(match(1:nTeams));
                teamRecord = zeros(length(teamList), 2);
                for teamIdx = 1:length(teamList)
                    match = find(d.league == uLeagues(whichLeague) & d.team == teamList(teamIdx));
                    teamRecord(teamIdx, 1) = d.wins(match(1));
                    teamRecord(teamIdx, 2) = d.games(match(2));
                end
                
                match = find(d.league ==  uLeagues(whichLeague));
                
                % figure window
                F = figure; clf; hold on;
                setFigure(F, [0.2 0.05 0.55 0.9], '');
                set(gcf, 'renderer', 'painters');
                
                set(gca, ...
                    'units'         ,  'normalized'           , ...
                    'position'      ,  [0.01 0.01 0.9 0.98]   ,...
                    'xlim'          ,  [0 nTeams+1]     , ...
                    'ylim'          ,  [-1 picksPerTeam+1]       , ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.01 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'clipping'      , 'off'                   , ...
                    'fontsize'      ,  fontSize                      );
                axis off;
                
                for teamIdx = 1:nTeams
                    team = find(teamList(teamIdx) == d.team(match(1:nTeams)));
                    str1 = sprintf('Team %d', teamIdx);
                    str2 = sprintf('%d-%d', teamRecord(team, 1), teamRecord(team, 2)-teamRecord(team, 1));
                    text(teamIdx, picksPerTeam+1, {str1, str2}, ...
                        'fontsize', fontSize+4, ...-
                        'hor', 'cen', ...
                        'vert', 'top');
                end
                
                rowIdx = picksPerTeam-1/2;
                direction = 1;
                for pickIdx = 1:length(match)
                    team = find(d.team(match(pickIdx)) == teamList);
                    %                 if ~isnan(team)
                    pos = find(sortIdx == teamList(team));
                    if ~isnan(d.role(match(pickIdx)))
                        rectangle('position', [pos-width/2 rowIdx-height/2 width height], ...
                            'curvature', [0.1 0.1], ...
                            'facecolor', [roleColorsFull{d.role(match(pickIdx))} 0.5], ...
                            'edgecolor', 'none');
                    end
                    str = d.playerNames{match(pickIdx)};
                    str = {str(1:strfind(str, ' ')-1), str(strfind(str, ' ')+1:end)};
                    text(pos, rowIdx, str, ...
                        'fontsize', fontSize, ...
                        'hor', 'cen', ...
                        'vert', 'mid');
                    switch direction
                        
                        case 1
                            
                            if pos > 1
                                [xaf,yaf] = ds2nfu(gca, [pos-1+width/2 pos-width/2], [rowIdx rowIdx]);
                                annotation('arrow', xaf, yaf, ...
                                    'headstyle', 'cback1', ...
                                    'headwidth', 5);
                            end
                            if pos == nTeams & pickIdx < length(match)
                                plot([pos+width/2 pos+width/2+gap], [rowIdx rowIdx], 'k-');
                                plot([pos+width/2+gap pos+width/2+gap], [rowIdx rowIdx-1], 'k-');
                                [xaf,yaf] = ds2nfu(gca, [pos+width/2+gap pos+width/2], [rowIdx rowIdx]-1);
                                annotation('arrow', xaf, yaf, ...
                                    'headstyle', 'cback1', ...
                                    'headwidth', 5);
                            end
                            
                        case -1
                            if pos > 1
                                [xaf,yaf] = ds2nfu(gca, [pos-width/2 pos-1+width/2],  [rowIdx rowIdx]);
                                annotation('arrow', xaf, yaf, ...
                                    'headstyle', 'cback1', ...
                                    'headwidth', 5);
                            elseif pos == 1 & pickIdx < length(match)
                                plot([pos-width/2 pos-width/2-gap],  [rowIdx rowIdx], 'k-');
                                plot([pos-width/2-gap pos-width/2-gap],  [rowIdx rowIdx-1], 'k-');
                                [xaf,yaf] = ds2nfu(gca, [pos-width/2-gap pos-width/2],  [rowIdx rowIdx]-1);
                                annotation('arrow', xaf, yaf, ...
                                    'headstyle', 'cback1', ...
                                    'headwidth', 5);
                                
                            end
                            %                     end
                    end
                    
                    if mod(pickIdx, nTeams) == 0
                        rowIdx = rowIdx - 1;
                        direction = -direction;
                    end
                    
                end
                
                for idx = 1:d.nRoles
                    H(idx) = patch([-10 -9 -9 -10], [-10 -10 -9 -9], 'k', ...
                        'facecolor', roleColorsFull{idx}, ...
                        'facealpha', 0.5, ...
                        'edgecolor', 'none');
                end
                L = legend(H, d.roleNames, ...
                    'box', 'off', ...
                    'fontsize', fontSize+8, ...
                    'location', 'east');
                set(L, 'pos', get(L, 'pos') + [0.1 0.05 0 0]);
                
                %             whichLeague
                %             pause
            end
            
        case 'handcuff'
            
            minLeagueSize = 6;
            maxLeagues = 1e3;
            requiredGames = 13;
            
            cuffPlayers = {'Adrian Peterson', 'Andre Ellington'; ...
                'Devonta Freeman', 'Tevin Coleman'; ...
                'Javorius Allen',	'Terrance West'; ...
                'LeSean McCoy',	 'Mike Tolbert'; ...
                'Christian McCaffrey',	 'Jonathan Stewart'; ...
                'Jordan Howard',	 'Tarik Cohen'; ...
                'Joe Mixon',	 'Jeremy Hill'; ...
                'Isaiah Crowell',	 'Duke Johnson'; ...
                'Ezekiel Elliott',	 'Darren McFadden'; ...
                'CJ Anderson',	 'Jamaal Charles'; ...
                'Ameer Abdullah',	 'Theo Riddick'; ...
                'Ty Montgomery',	 'Aaron Jones'; ...
                'Lamar Miller',	 'D''Onta Foreman'; ...
                'Frank Gore',	 'Marlon Mack'; ...
                'Leonard Fournette',	 'Chris Ivory'; ...
                'Kareem Hunt',	 'Charcandrick West'; ...
                'Melvin Gordon',	 'Austin Ekeler'; ...
                'Todd Gurley',	 'Malcolm Brown'; ...
                'Jay Ajayi',	 'Kenyan Drake'; ...
                'Jerick McKinnon',	 'Latavius Murray'; ...
                'Mike Gillislee',	 'James White'; ...
                'Mark Ingram',	 'Alvin Kamara'; ...
                'Orleans Darkwa',	 'Wayne Gallman'; ...
                'Bilal Powell',	 'Matt Forte'; ...
                'Marshawn Lynch',	 'Deandre Washington'; ...
                'LeGarrette Blount',	 'Wendell Smallwood'; ...
                'Le''Veon Bell',	 'James Conner'; ...
                'Thomas Rawls',	 'Eddie Lacy'; ...
                'Carlos Hyde',	 'Matt Breida'; ...
                'Doug Martin',	 'Jacquizz Rodgers'; ...
                'DeMarco Murray',	 'Derrick Henry'; ...
                'Rob Kelley',	 'Chris Thompson'};
            
            [nCuff, ~] = size(cuffPlayers);
            
            % 1 = player1 in league, 2 = player2 in league, 3 = both in league,
            % 4 = both in same team, 5= chance
            m = zeros(nCuff, 5);
            nLeagues = 0; leagueSize = nan(maxLeagues, 1);
            winsCuff = []; gamesCuff = [];
            winsNoCuff = []; gamesNoCuff = [];
            leagueList = unique(d.league);
            for leagueIdx = 1:length(leagueList)
                match = find(d.league == leagueList(leagueIdx));
                
                go = true;
                tmp = length(unique(d.team(match)));
                if go
                    nLeagues = nLeagues + 1;
                    leagueSize(nLeagues)= tmp;
                    for idx = 1:nCuff
                        match1 = find(strcmp(cuffPlayers{idx, 1}, d.playerNames(match)));
                        if ~isempty(match1)
                            match1 = match1(find(d.auto(match1) ~= 1));
                        end
                        if ~isempty(match1)
                            m(idx, 1) = m(idx, 1) + 1;
                            player1 = d.team(match(match1(1)));
                        end
                        match2 = find(strcmp(cuffPlayers{idx, 2}, d.playerNames(match)));
                        if ~isempty(match2)
                            match2 = match2(find(d.auto(match2) ~= 1));
                        end
                        if ~isempty(match2)
                            m(idx, 2) = m(idx, 2) + 1;
                            player2 = d.team(match(match2(1)));
                        end
                        if ~isempty(match1) & ~isempty(match2)
                            m(idx, 3) = m(idx, 3) + 1;
                            m(idx, 5) = m(idx, 5) + 1/leagueSize(nLeagues);
                            if player1 == player2
                                m(idx, 4) = m(idx, 4) + 1;
                                winsCuff = [winsCuff d.wins(match(match1(1)))];
                                gamesCuff = [gamesCuff d.games(match(match1(1)))];
                            else
                                winsNoCuff = [winsNoCuff d.wins(match(match1(1)))];
                                gamesNoCuff = [gamesNoCuff d.games(match(match1(1)))];
                            end
                        end
                    end
                    
                end
            end
            leagueSize = leagueSize(1:nLeagues);
            keepCuff = find(ismember(gamesCuff, [13 14]));
            winsCuff = winsCuff(keepCuff);
            gamesCuff = gamesCuff(keepCuff);
            keepNoCuff = find(ismember(gamesNoCuff, [13 14]));
            winsNoCuff = winsNoCuff(keepNoCuff);
            gamesNoCuff = gamesNoCuff(keepNoCuff);
            
            % dump latex table
            fid = fopen('handcuffTable.txt', 'w');
            fprintf(fid, '\\begin{table}\n');
            fprintf(fid, '\\begin{center}\n');
            % fprintf(fid, '\\resizebox{\\textwidth}{!}{\n');
            fprintf(fid, '\\begin{tabular}{llrrrrr}\n');
            fprintf(fid, '\\toprule\n');
            fprintf(fid, ' First Player & Second Player & First & Second & Both & Same & Chance \\\\ \n');
            fprintf(fid, '\\hline\n');
            for idx = 1:nCuff
                fprintf(fid, sprintf('%s & %s & %d & %d & %d & %d & %d \\\\\\\\ \n', cuffPlayers{idx, :}, round(m(idx, :))));
            end
            fprintf(fid, '\\bottomrule\n');
            fprintf(fid, '\\end{tabular}\n');
            fprintf(fid, '\\end{center}\n');
            fprintf(fid, '\\end{table}\n');
            
            % test mean win proportion difference
            % which engine to use
            engine = 'jags';
            preLoad = true;
            
            dataName = 'handcuff2017';
            
            % graphical model script
            modelName = 'hierarchicalRateDifference_1';
            
            % parameters to monitor
            params = {'mu', 'sigma', 'delta', 'deltaPrior'};
            
            % MCMC properties
            nChains    = 8;     % nuber of MCMC chains
            nBurnin    = 1e3;   % number of discarded burn-in samples
            nSamples   = 1e3;   % number of collected samples
            nThin      = 1;     % number of samples between those collected
            doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains
            
            % assign MATLAB variables to the observed nodes
            data = struct(...
                'y1' , winsCuff , ...
                'n1' , gamesCuff, ...
                'y2' , winsNoCuff , ...
                'n2' , gamesNoCuff, ...
                'p1' , length(winsCuff), ...
                'p2' , length(winsNoCuff)    );
            
            % generator for initialization
            generator = @()struct('muGrand', rand, 'sigma', [0.1 0.1]);
            
            % Sample using Trinity
            fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
            
            if preLoad && isfile(sprintf('storage/%s', fileName))
                fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
                load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
            else
                tic; % start clock
                [stats, chains, diagnostics, info] = callbayes(engine, ...
                    'model'           , sprintf('%s_%s.txt', modelName, engine)   , ...
                    'data'            , data                                      , ...
                    'outputname'      , 'samples'                                 , ...
                    'init'            , generator                                 , ...
                    'datafilename'    , modelName                                 , ...
                    'initfilename'    , modelName                                 , ...
                    'scriptfilename'  , modelName                                 , ...
                    'logfilename'     , sprintf('tmp/%s', modelName)              , ...
                    'nchains'         , nChains                                   , ...
                    'nburnin'         , nBurnin                                   , ...
                    'nsamples'        , nSamples                                  , ...
                    'monitorparams'   , params                                    , ...
                    'thin'            , nThin                                     , ...
                    'workingdir'      , sprintf('tmp/%s', modelName)              , ...
                    'verbosity'       , 0                                         , ...
                    'saveoutput'      , true                                      , ...
                    'allowunderscores', true                                      , ...
                    'parallel'        , doParallel                                );
                fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing
                fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
                if ~isfolder('storage')
                    !mkdir storage
                end
                save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
                
                % convergence of each parameter
                disp('Convergence statistics:')
                grtable(chains, 1.05)
                
                % basic descriptive statistics
                disp('Descriptive statistics for all chains:')
                codatable(chains);
                
            end
            
            % convergence of each parameter
            disp('Convergence statistics:')
            grtable(chains, 1.05)
            
            % basic descriptive statistics
            disp('Descriptive statistics for all chains:')
            codatable(chains);
            
            % Savage-Dickey on delta
            
            % constants
            critical = 0;
            fontSize = 18;
            lineWidth = 2;
            markerSize = 10;
            posteriorColor = pantone.ClassicBlue;
            priorColor = pantone.Titanium;
            
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.6 0.4], '');
            
            for idx = 1:2
                
                switch idx
                    case 1
                        binHi = 1; binLo = -1; binWidth = 0.005; binTick = 0.5;
                    case 2
                        binLo = -0.01; binHi = 0.01; binWidth = 0.001; binTick = 0.01;
                end
                binCenter = binLo:binWidth:binHi;
                [~, criticalIdx] = min(abs(binCenter - critical));
                
                subplot(1, 2, idx); hold on;
                
                % axis
                set(gca, ...
                    'ycolor'     , 'none'                , ...
                    'xlim'       , [binLo binHi]         , ...
                    'xtick'      , binLo:binTick:binHi   , ...
                    'box'        , 'off'                 , ...
                    'tickdir'    , 'out'                 , ...
                    'layer'      , 'top'                 , ...
                    'ticklength' , [0.02 0]              , ...
                    'layer'      , 'top'                 , ...
                    'fontsize'   , fontSize              );
                set(gca, 'pos', get(gca, 'pos') + [-0.05 0.075 0 0]);
                
                % labels
                xlabel('Difference in Rates', 'fontsize', fontSize+2);
                
                densityPrior = histcounts(chains.deltaPrior(:), ...
                    'binlimits'     ,[binLo-binWidth/2 binHi+binWidth/2]    , ...
                    'binwidth'      , binWidth         , ...
                    'normalization' , 'pdf'            );
                plot(binCenter, densityPrior, '--', ...
                    'color'     , priorColor  , ...
                    'linewidth' , lineWidth   );
                
                densityPosterior = histcounts(chains.delta(:), ...
                    'binlimits'     , [binLo-binWidth/2 binHi+binWidth/2]    , ...
                    'binwidth'      , binWidth         , ...
                    'normalization' , 'pdf'            );
                plot(binCenter, densityPosterior, '-', ...
                    'color'     , posteriorColor  , ...
                    'linewidth' , lineWidth   );
                
                
                if idx == 1
                    legend('prior', 'posterior', ...
                        'box', 'off', ...
                        'fontsize', fontSize, ...
                        'location', 'northwest', ...
                        'autoupdate', 'off');
                end
                
                BF = exp(log(densityPrior(criticalIdx)) - log(densityPosterior(criticalIdx)));
                if idx == 2
                    if BF > 1
                        str = sprintf(' BF_{10} = %1.1f', BF);
                    else
                        str = sprintf(' BF_{01} = %1.1f', 1/BF);
                        
                    end
                    plot(binCenter(criticalIdx)*ones(1, 2), [densityPrior(criticalIdx) densityPosterior(criticalIdx)], 'k-');
                    text(binCenter(criticalIdx), mean([densityPrior(criticalIdx) densityPosterior(criticalIdx)]), str, ...
                        'fontsize', fontSize, ...
                        'horizontal', 'left');
                end
                
                plot(binCenter(criticalIdx), densityPrior(criticalIdx), 'o', ...
                    'linewidth'   , lineWidth, ...
                    'markeredgecolor',  priorColor, ...
                    'markerfacecolor', 'w', ...
                    'markersize', markerSize);
                plot(binCenter(criticalIdx), densityPosterior(criticalIdx), 'o', ...
                    'linewidth'   , lineWidth, ...
                    'markeredgecolor',  posteriorColor, ...
                    'markerfacecolor', 'w', ...
                    'markersize', markerSize);
                
            end
            % print
            if doPrint
                warning off;
                print(sprintf('figures/%s_%s_SavageDickey.png', modelName, dataName), '-dpng');
                print(sprintf('figures/%s_%s_SavageDickey.pdf', modelName, dataName), '-dpdf');
                print(sprintf('figures/%s_%s_SavageDickey.eps', modelName, dataName), '-depsc');
                warning on;
            end
            
        case 'playersInset'
            
            maxPlayersShow = 50;
            maxPlayers = 150;
            maxSlots = 120;
            minPicks = 50;
            fontSize = 12;
            bounds = [25 75];
            winThreshold = [0.465 0.57];
            
            m = zeros(maxPlayers, maxSlots);
            winProp = cell(maxPlayers, 1);
            pickMean = zeros(maxPlayers, 1);
            for idx = 0:d.nPlayers-1
                match = find(d.player == idx & d.team ~= -1 & d.auto ~= 1);
                count = histcounts(d.pick(match), 0.5:maxSlots+1.5);
                m(idx+1, 1:maxSlots) = count(1:end-1);
                if sum(count) >= minPicks
                    pickMean(idx+1) = mean(d.pick(match));
                else
                    pickMean(idx+1) = inf;
                end
                winProp{idx+1} = d.wins(match)./d.games(match);
            end
            
            [pickMean, sortIdx] = sort(pickMean, 'ascend');
            
            yTickLabel = cell(maxPlayersShow, 1);
            role = nan(maxPlayersShow, 1);
            for idx = 1:maxPlayersShow
                match = find(d.player == sortIdx(idx)-1);
                yTickLabel{idx} = d.playerNames{match(1)};
                role(idx) = d.role(match(1));
            end
            m = m(sortIdx, :);
            winProp = winProp(sortIdx);
            
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.1 0.4 0.8], '');
            
            set(gca, ...
                'units'         ,  'normalized'           , ...
                'position'      ,  [0.245 0.1 0.735 0.85]     ,...
                'xlim'          ,  [1 maxSlots]         , ...
                'xtick'         ,  [1 20:20:maxSlots]     , ...
                'ylim'          ,  [1 maxPlayersShow]       , ...
                'ydir'          ,  'reverse'              , ...
                'ytick'         ,  1:maxPlayersShow           , ...
                'yticklabel'    ,  yTickLabel   , ...
                'box'           ,  'off'                  , ...
                'tickdir'       ,  'out'                  , ...
                'ticklength'    ,  [0.0075 0]               , ...
                'layer'         ,  'top'                  , ...
                'fontsize'      ,  fontSize                      );
            xlabel('Draft Position', 'fontsize', fontSize+5);
            Raxes(gca, 0.01, 0.02);
            
            % plot
            scale = 25/max(sqrt(m(:)));
            for playerIdx = 1:maxPlayersShow
                for slotIdx = 1:maxSlots
                    if m(playerIdx, slotIdx) > 0
                        plot(slotIdx, playerIdx, 'ko', ...
                            'markersize', scale*sqrt(m(playerIdx, slotIdx))            , ...
                            'markerfacecolor' , roleColorsFull{role(playerIdx)} , ...
                            'markeredgecolor' , 'none'                                          );
                    end
                end
            end
            
            meanWins = cellfun(@mean, winProp);
            seWins = cellfun(@std, winProp)./sqrt(cellfun(@length, winProp));

            matchLag = find(pickMean~=Inf & ~isnan(pickMean) & ~isnan(meanWins));
            AX = axes; hold on;
            set(AX, ...
                'units'         ,  'normalized'           , ...
                'position'      ,  [0.625 0.625 0.35 0.35]     ,...
                'xlim'          ,  [1 maxPlayers]         , ...
                'xtick'         ,  [1 maxPlayers]   , ...
                'ylim'          ,  [0.45 0.61]       , ...
                'ytick'         ,  [0.45:0.05:0.6]          , ...
                'box'           ,  'off'                  , ...
                'tickdir'       ,  'out'                  , ...
                'ticklength'    ,  [0.02 0]               , ...
                'layer'         ,  'top'                  , ...
                'fontsize'      ,  fontSize                      );
            xlabel('Pick Position', 'fontsize', fontSize+3);
            ylabel('Winning Proportion', 'fontsize', fontSize+3);
            Raxes(gca, 0.02, 0.02);
            
            for idx = 1:length(matchLag)
                match = find(d.player == sortIdx(matchLag(idx))-1);
                if ~isnan(d.role(match(1)))
                    plot(pickMean(matchLag(idx)), meanWins(matchLag(idx)), 'o', ...
                        'markersize', 6, ...
                        'markerfacecolor', roleColorsFull{d.role(match(1))}, ...
                        'markeredgecolor', 'w');
                    if meanWins(matchLag(idx)) < winThreshold(1) | meanWins(matchLag(idx)) > winThreshold(2)
                        str = split(d.playerNames{match(1)}, ' ');
                        text(pickMean(matchLag(idx)), meanWins(matchLag(idx)), [' ' str{2}], ...
                            'fontsize', fontSize);
                    end
                end
            end            
            
            
            
        case 'copycat3'
            
            % plot probability a role is taken in each draft position
            % unconditionally, and conditional on same role just selected
            
            maxPick =  120; %max(d.pick);
            fontSize = 16;
            threshold = 20; % how many before taking fraction
            
            if exist('copycat3mSave.mat')
                load copycat3mSave m
            else
                m = zeros(maxPick, d.nRoles, 6);
                for roleIdx = 1:d.nRoles
                    for pickIdx = 1:maxPick
                        match = find(d.pick == pickIdx-1  & d.auto == 0 & d.team ~= -1);
                        m(pickIdx, roleIdx, 2) = m(pickIdx, roleIdx, 2) + length(match);
                        match = find(d.pick == pickIdx-1 & d.role == roleIdx & d.auto == 0 & d.team ~= -1);
                        m(pickIdx, roleIdx, 1) = m(pickIdx, roleIdx, 1) + length(match);
                        if pickIdx > 1
                            match = find(d.pick == pickIdx-1 & d.role == roleIdx  & d.auto == 0 & d.team ~= -1 & d.team ~= circshift(d.team, -1));
                            m(pickIdx, roleIdx, 4) = m(pickIdx, roleIdx, 4) + length(match);
                            match = find(d.pick == pickIdx-1 & d.role == roleIdx & circshift(d.role, -1) == roleIdx  & d.auto == 0 & d.team ~= -1 & d.team ~= circshift(d.team, -1));
                            m(pickIdx, roleIdx, 3) = m(pickIdx, roleIdx, 3) + length(match);
                        end
                        match = find(d.pick == pickIdx-1  & d.auto == 0 & d.team ~= -1 & d.role == roleIdx);
                        for idx = 1:length(match)
                            currentTeam = d.team(match(idx));
                            currentLeague = d.league(match(idx));
                            matchTeam = find(d.team == currentTeam & d.league == currentLeague & d.pick < pickIdx-1);
                            if ~isempty(matchTeam)
                                m(pickIdx, roleIdx, 6) = m(pickIdx, roleIdx, 6) + 1;
                                if d.role(matchTeam(end)) == roleIdx
                                    % [pickIdx roleIdx currentTeam currentLeague]
                                    m(pickIdx, roleIdx, 5) = m(pickIdx, roleIdx, 5) + 1;
                                end
                            end
                        end
                        [roleIdx pickIdx]
                    end
                end
                save copycat3mSave m
            end
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.1 0.5 0.8], '');
            set(gcf, 'renderer', 'painters');
            
            for roleIdx = 1:d.nRoles
                
                subplot(d.nRoles/2, 2, roleIdx); hold on;
                % composition axes
                set(gca, ...
                    'xlim'          ,  [1 maxPick]     , ...
                    'xtick'         ,  [1 20:20:maxPick]            , ...
                    'ylim'          ,  [0 1]       , ...
                    'ytick'         , 0:0.2:1          , ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.01 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'XTickLabelRotation', 0 , ...
                    'fontsize'      ,  fontSize                      );
                
                val = (m(:, roleIdx, 1)+1)./(m(:, roleIdx, 2) + 2);
                patch([1:maxPick maxPick 1], [val' 0 0], 'k',  ...
                    'facecolor', roleColorsFull{roleIdx}, ...
                    'edgecolor', roleColorsFull{roleIdx}, ...
                    'facealpha', 0.5);
                
                
                val = (m(:, roleIdx, 3)+1)./(m(:, roleIdx, 4) + 2);
                keep = find(m(:, roleIdx, 4) > threshold);
                H(1)=  plot(keep, val(keep), '-', ...
                    'markerfacecolor', 'w', ...
                    'markersize', 3, ...
                    'markeredgecolor', pantone.ClassicBlue, ...
                    'color', pantone.ClassicBlue, ...
                    'linewidth', 2);
                
                text(1, 1, d.roleNames{roleIdx}, ...
                    'vert', 'top', ...
                    'hor', 'left', ...
                    'fontsize', fontSize+2);
                
                Raxes(gca, 0.005, 0.005);
                
            end
            
            [~, H(1)] = suplabel('Draft Position');
            [~, H(2)] = suplabel('Choice Probability', 'y');
            set(H, 'fontsize', fontSize + 4);
            
        case 'positionChoice'
            
            maxPick =  150;
            fontSize = 16;
            threshold = 20; % how many before taking fraction
            
            m = zeros(maxPick, d.nRoles, 2);
            for roleIdx = 1:d.nRoles
                for pickIdx = 1:maxPick
                    match = find(d.pick == pickIdx-1  & d.auto == 0 & d.team ~= -1);
                    m(pickIdx, roleIdx, 2) = m(pickIdx, roleIdx, 2) + length(match);
                    match = find(d.pick == pickIdx-1 & d.role == roleIdx & d.auto == 0 & d.team ~= -1);
                    m(pickIdx, roleIdx, 1) = m(pickIdx, roleIdx, 1) + length(match);
                end
            end
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.1 0.6 0.6], '');
            set(gcf, 'renderer', 'painters');
            
            % composition axes
            set(gca, ...
                'xlim'          ,  [1 maxPick]     , ...
                'xtick'         ,  [1 30:30:maxPick]            , ...
                'ylim'          ,  [0 1]       , ...
                'ytick'         , 0:0.2:1          , ...
                'box'           ,  'off'                  , ...
                'tickdir'       ,  'out'                  , ...
                'ticklength'    ,  [0.01 0]               , ...
                'layer'         ,  'top'                  , ...
                'fontsize'      ,  fontSize                      );
            xlabel('Draft Position', 'fontsize', fontSize+4);
            ylabel('Choice Probability', 'fontsize', fontSize+4);
            
            
            val = (m(:, :, 1))./(m(:, :, 2));
            val = val./repmat(sum(val, 2), [1 d.nRoles]);
            [~, sortIdx] = sort(sum(val), 'descend');
            H = bar(1:maxPick, val(:, sortIdx), 1, 'stacked');
            for i = 1:numel(H)
                set(H(i), 'facecolor', roleColorsFull{sortIdx(i)}, ...
                    'facealpha', 0.7, ...
                    'edgecolor', pantone.Titanium, ...
                    'linewidth', 0.5);
                
                for i = 1:d.nRoles
                    if i == 1
                        text(maxPick, sum(val(maxPick, sortIdx(1:i)))/2, [' ' d.roleNames{sortIdx(i)}], ...
                            'vert', 'mid', ...
                            'hor', 'left', ...
                            'fontsize', fontSize);
                    else
                        text(maxPick, sum(val(maxPick, sortIdx(1:i-1))) + val(maxPick, sortIdx(i))/2, [' ' d.roleNames{sortIdx(i)}], ...
                            'vert', 'mid', ...
                            'hor', 'left', ...
                            'fontsize', fontSize);
                    end
                end
            end
            
            Raxes(gca, 0.005, 0.005);
            
        case 'cumulativePicks'
            
            % panels with distribution how many of each role have been chosen by a team
            % after 1, 2, ... picks
            
            maxPick = 15;
            maxCount = 20;
            m = zeros(d.nRoles, maxCount, maxPick);
            scale = 1;
            maxHowMany = 0;
            fontSize = 14;
            
            % from https://www.nfl.com/news/the-best-fantasy-football-draft-strategy-for-2017-0ap3000000815105
            % set upper and lower bounds on each role after each pick
            bounds{1} = [0 0; 0 0; 0 0; 0 0; 0 1; 0 1; 0 1; 0 2; 0 2; 1 3; 1 3; 1 4; 1 4; 1 4; 1 4]; % QB
            bounds{2} = [0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 7; 0 8; 0 9; 0 10; 0 11; 0 13; 0 13; 0 13]; % RB
            bounds{3} = [0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 0 7; 0 7; 0 8; 0 9; 0 10; 0 11; 0 13; 0 13; 0 13]; % WR
            bounds{4} = [0 0; 0 0; 0 0; 0 0; 0 1; 0 1; 0 2; 0 3; 0 3; 1 4; 1 4; 1 5; 1 5; 1 5; 1 5]; % TE
            bounds{5} = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 1; 1 1]; % K
            bounds{6} = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 1; 1 1]; % DS
            
            if exist('cumulativePicks.mat')
                load cumulativePicks m maxHowMany
            else
                for leagueIdx = 1:d.nLeagues
                    match = find(d.league == leagueIdx-1);
                    teamList = unique(d.team(match));
                    teamList = teamList(find(teamList >= 0));
                    for teamIdx = 1:length(teamList)
                        match = find(d.league == leagueIdx-1 & d.team == teamList(teamIdx));
                        match = match(1:min(end, maxPick));
                        for idx = 1:length(match)
                            for roleIdx = 1:d.nRoles
                                howMany = sum(d.role(match(1:idx)) == roleIdx);
                                maxHowMany = max(maxHowMany, howMany);
                                m(roleIdx,  howMany+1, idx) = m(roleIdx,  howMany+1, idx) + 1;
                            end
                        end
                    end
                end
                m = m(:, 1:maxHowMany+1, :);
                save cumulativePicks m maxHowMany
            end
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.375 0.6], '');
            
            for roleIdx = 1:d.nRoles
                
                subplot(d.nRoles/2, 2, roleIdx); hold on;
                % composition axes
                set(gca, ...
                    'xlim'          ,  [0 maxPick+1]     , ...
                    'xtick'         ,  [1 5 10 15]            , ...
                    'ylim'          ,  [-1 8]       , ...
                    'ytick'         , 0:1:7          , ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.02 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'fontsize'      ,  fontSize                      );
                set(gca, 'pos', get(gca, 'pos') + [0 0.05 0 0]);
                
                for pickIdx = 1:maxPick
                    for hmIdx = 0:maxHowMany
                        val = sqrt(m(roleIdx, hmIdx+1, pickIdx)/sum(m(roleIdx, :, pickIdx)));
                        if val > 0
                            rectangle('position', [pickIdx-scale*val/2 hmIdx-scale*val/2 scale*val scale*val], ...
                                'facecolor', [roleColorsFull{roleIdx} 1], ...
                                'edgecolor', 'w', ...
                                'curvature', [1 1 ]);
                            
                        end
                    end
                    tmp(1) = plot([pickIdx-scale/2 pickIdx+scale/2], bounds{roleIdx}(pickIdx, 1) * ones(1, 2)-scale/2, '-');
                    tmp(2) =plot([pickIdx-scale/2 pickIdx+scale/2], bounds{roleIdx}(pickIdx, 2) * ones(1, 2)+scale/2, '-');
                    if pickIdx > 1
                        tmp(3) =plot(pickIdx-1+scale/2*ones(1,2), bounds{roleIdx}(pickIdx-1:pickIdx, 1)-scale/2, '-');
                        tmp(4) =plot(pickIdx-1+scale/2*ones(1,2), bounds{roleIdx}(pickIdx-1:pickIdx, 2)+scale/2, '-');
                    end
                    set(tmp, 'color', pantone.Titanium);
                    
                end
                
                text(1, max(get(gca, 'ylim')), d.roleNames{roleIdx}, ...
                    'vert', 'top', ...
                    'hor', 'left', ...
                    'fontsize', fontSize+2);
                
                Raxes(gca, 0.005, 0.005);
                
            end
            
            [~, H(1)] = suplabel('Draft Position');
            [~, H(2)] = suplabel('Cumulative Picks', 'y');
            set(H, 'fontsize', fontSize + 4);
            
        case 'descriptiveStats'
            
            nHumanTeams = 0;
            nTeams = 0;
            leagueList = unique(d.league);
            m = zeros(20, 1);
            for leagueIdx = 1:length(leagueList)
                match = find(d.league == leagueList(leagueIdx));
                nHumanTeams = nHumanTeams + d.leagueHumans(match(1));
                nTeams = nTeams + d.leagueTeams(match(1));
                diff = d.leagueTeams(match(1)) - d.leagueHumans(match(1));
                m(diff+1) = m(diff+1) + 1;
            end
            fprintf('There are %d leagues, with a total of %d teams, of which %d are human players and %d are bots.\n', ...
                d.nLeagues, nTeams, nHumanTeams, nTeams - nHumanTeams);
            
            match = find(d.team >= 0);
            t = tabulate(d.auto(match));
            fprintf('There are %d picks made by teams owned by people. %1.0f%% are autopicks.\n', length(match), t(2,3));
            
            fprintf('%1.0f%% of leagues have all human teams, %1.0f%% have one bot, %1.0f%% have two bots, and %1.0f%% have three bots.\n', ...
                m(1:4)./sum(m)*100);
            
            
        case 'fandom2'
            
            
            maxTeams = 2e4;
            sizeList = 14:20; % sizeList = 15;
            nRows = 8; nCols = 4;
            fontSize = 16;
            
            t = tabulate(d.nflTeam);
            teamBase = t(:, 3)/100;
            [~, baseIdx] = sort(teamBase, 'descend');
            
            nflCounts = zeros(length(sizeList), numel(d.nflTeamNames), max(sizeList)+1);
            
            leagueList = unique(d.league);
            for  leagueIdx= 1:length(leagueList)
                match = find(d.league == leagueList(leagueIdx));
                teamList = unique(d.team(match));
                teamList = teamList(find(teamList >= 0));
                for teamIdx = 1:length(teamList)
                    match = find(d.league == leagueList(leagueIdx) & d.team == teamList(teamIdx));
                    sizeMatch = find(length(match) == sizeList);
                    counts = histcounts(d.nflTeam(match), 0.5:numel(d.nflTeamNames)+0.5);
                    for idx = 1:length(counts)
                        nflCounts(sizeMatch, idx, counts(idx)+1) = nflCounts(sizeMatch, idx, counts(idx)+1) + 1;
                    end
                end
            end
            
            for sizeIdx = 1:length(sizeList)
                
                size = sizeList(sizeIdx);
                
                % figure window
                F = figure; clf; hold on;
                setFigure(F, [0.2 0.1 0.4 0.7], '');
                set(gcf, 'Renderer', 'painters');
                
                for nflIdx = 1:numel(d.nflTeamNames)
                    
                    subplot(nRows, nCols, nflIdx); hold on;
                    % composition axes
                    set(gca, ...
                        'xlim'          ,  [-0.5 8.5]     , ...
                        'xtick'         ,  0:1:8          , ...
                        'ylim'          ,  [0 1]       , ...
                        'ytick'         , 0:0.2:1        , ...
                        'xticklabel'    , {'0', '', '', '', '', '', '', '', '8'}, ...
                        'yticklabel'    , {'0', '', '', '', '', '1'}, ...
                        'XTickLabelRotation', 0, ...
                        'box'           ,  'off'                  , ...
                        'tickdir'       ,  'out'                  , ...
                        'ticklength'    ,  [0.03 0]               , ...
                        'layer'         ,  'top'                  , ...
                        'fontsize'      ,  fontSize                      );
                    if nflIdx ~= (numel(d.nflTeamNames) - nCols + 1)
                        set(gca, ...
                            'xticklabel'    , [], ...
                            'yticklabel' , []);
                    end
                    set(gca, 'pos', get(gca, 'pos').*[1 1 1 1] + [0 0.05 0 0]);
                    Raxes(gca, 0.01, 0.005);
                    
                    hi = zeros(size+1, 1);
                    theta = teamBase(baseIdx(nflIdx));
                    for idx = 0:size
                        hi(idx+1) = nchoosek(size, idx) * theta^idx * (1-theta)^(size-idx);
                    end
                    hi = hi/sum(hi);
                    
                    H(1) = bar(0:size, hi, 1, 'k');
                    set(H(1), 'facecolor', pantone.Sandstone, ...
                        'edgecolor',  pantone.Sandstone, ...
                        'facealpha', 0.5);
                    
                    count = squeeze(nflCounts(sizeIdx, baseIdx(nflIdx), :));
                    count = count/sum(count);
                    keep = find(count > 0);
                    H(2) = bar(keep-1, count(keep), 0.4, 'k');
                    set(H(2), 'facecolor', pantone.ClassicBlue, ...
                        'edgecolor', 'none');
                    
                    text(8, 1, lower(d.nflTeamNames{baseIdx(nflIdx)}), ...
                        'fontsize', fontSize, ...
                        'vert', 'top', ...
                        'hor', 'right');
                    
                end
                
                [~, H(1)] = suplabel('Number of Players in Team');
                [~, H(2)] = suplabel('Probability', 'y');
                set(H(2), 'pos', get(H(2), 'pos') + [0.035 0 0]);
                set(H(1), 'pos', get(H(1), 'pos') + [0 0.02 0]);
                set(H, 'fontsize', fontSize + 2);
                
                if doPrint
                    print(sprintf('figures/%s_%d.png', aList{aListi}, size), '-dpng', '-r600');
                    print(sprintf('%s/%s_%d.eps', figDir, aList{aListi}, size), '-depsc');
                end
                
            end
            return
            
        case 'positionPerformance'
            
            teamList = 6:2:12;
            fontSize = 18;
            scale = 0.02; gap = 0.05;
            
            wins = cell(length(teamList), max(teamList));
            
            leagueList = unique(d.league);
            for teamIdx = 1:length(teamList)
                for  leagueIdx= 1:length(leagueList)
                    match = find(d.league == leagueList(leagueIdx) & d.leagueTeams == teamList(teamIdx));
                    if ~isempty(match)
                        for posIdx = 1: teamList(teamIdx)
                            wins{teamIdx, posIdx} = [wins{teamIdx, posIdx} d.wins(match(posIdx))/d.games(match(posIdx))];
                        end
                    end
                end
            end
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.6 0.3], '');
            
            left = 0.1;
            
            for teamIdx = 1:length(teamList)
                subplot(1, length(teamList), teamIdx); hold on;
                % composition axes
                set(gca, ...
                    'units'         , 'normalized', ...
                    'pos'           , [left 0.225 teamList(teamIdx)*scale 0.7], ...
                    'xlim'          ,  [1 teamList(teamIdx)]     , ...
                    'xtick'         ,  1:teamList(teamIdx)        , ...
                    'ylim'          ,  [0.3 0.7]       , ...
                    'ytick'         , 0.3:0.1:0.7        , ...
                    'XTickLabelRotation', 0, ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.02 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'fontsize'      ,  fontSize                      );
                left = left + teamList(teamIdx)*scale + gap;
                if teamIdx ~= 1
                    set(gca, 'yticklabel', []);
                end
                text(teamList(teamIdx), max(get(gca, 'ylim')), {sprintf('%d team', teamList(teamIdx)), 'leagues'}, ...
                    'fontsize', fontSize, ...
                    'vert', 'top', ...
                    'hor', 'right');
                Raxes(gca, 0.02, 0.01);
                
                meanWins = cellfun(@nanmean, wins(teamIdx,  1:teamList(teamIdx) ));
                seWins = cellfun(@nanstd, wins(teamIdx,  1:teamList(teamIdx) ))./sqrt(cellfun(@length, wins(teamIdx,  1:teamList(teamIdx) )));
                errorbar(1:teamList(teamIdx), meanWins, seWins, ...
                    'marker', 'o', ...
                    'markersize', 6, ...
                    'markerfacecolor', winsColor, ...
                    'markeredgecolor', 'w', ...
                    'color', winsColor);
            end
            
            [~, H(1)] = suplabel('Drafting Position in League');
            [~, H(2)] = suplabel('Win Proportion', 'y');
            set(H(2), 'pos', get(H(2), 'pos') + [0.025 0 0]);
            set(H, 'fontsize', fontSize + 4);
            
        case 'copycatConsequences'
            
            % look at whether made above base rate copycat choices
            % against win percentage
            
            % first of role, anywhere in draft
            roleList  = [1 2 3 4 5 6];
            for roleIdx = 1:length(roleList)
                role = roleList(roleIdx);
                winsCopy = []; winsNoCopy = [];
                gamesCopy = []; gamesNoCopy = [];
                leagueList = unique(d.league);
                for  leagueIdx= 1:length(leagueList)
                    match = find(d.league == leagueList(leagueIdx));
                    teamList = unique(d.team(match));
                    teamList = teamList(find(teamList >= 0));
                    for teamIdx = 1:length(teamList)
                        match = find(d.league == leagueList(leagueIdx) & d.team == teamList(teamIdx));
                        spot = find(d.role(match) == role);
                        if ~isempty(spot)
                            if d.pick(match(spot(1))) > 1 & d.auto(match(spot(1))) == 0
                                if d.role(match(spot(1))-1) == role
                                    winsCopy = [winsCopy d.wins(match(1))];
                                    gamesCopy = [gamesCopy d.games(match(1))];
                                else
                                    winsNoCopy = [winsNoCopy d.wins(match(1))];
                                    gamesNoCopy = [gamesNoCopy d.games(match(1))];
                                end
                            end
                        end
                    end
                end
                fprintf('For first pick of role %s\ncopying teams mean win is %1.4f\nnon-copying teams is %1.4f\n---\n', ...
                    d.roleNames{role}, nanmean(winsCopy./gamesCopy), nanmean(winsNoCopy./gamesNoCopy));
                
                keepCopy = find(ismember(gamesCopy, [13 14]));
                winsCopy = winsCopy(keepCopy);
                gamesCopy = gamesCopy(keepCopy);
                keepNoCopy = find(ismember(gamesNoCopy, [13 14]));
                winsNoCopy = winsNoCopy(keepNoCopy);
                gamesNoCopy = gamesNoCopy(keepNoCopy);
                % test mean win proportion difference
                % which engine to use
                engine = 'jags';
                preLoad = true;
                
                dataName = 'copycat2017';
                
                % graphical model script
                modelName = 'hierarchicalRateDifference_1';
                
                % parameters to monitor
                params = {'mu', 'sigma', 'delta', 'deltaPrior'};
                
                % MCMC properties
                nChains    = 8;     % nuber of MCMC chains
                nBurnin    = 1e3;   % number of discarded burn-in samples
                nSamples   = 1e3;   % number of collected samples
                nThin      = 1;     % number of samples between those collected
                doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains
                
                % assign MATLAB variables to the observed nodes
                data = struct(...
                    'y1' , winsCopy , ...
                    'n1' , gamesCopy, ...
                    'y2' , winsNoCopy , ...
                    'n2' , gamesNoCopy, ...
                    'p1' , length(winsCopy), ...
                    'p2' , length(winsNoCopy)    );
                
                % generator for initialization
                generator = @()struct('muGrand', rand, 'sigma', [0.1 0.1]);
                
                % Sample using Trinity
                fileName = sprintf('%s_%s_%s_%s.mat', modelName, dataName, d.roleNames{roleIdx}, engine);
                
                if preLoad && isfile(sprintf('storage/%s', fileName))
                    fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
                    load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
                else
                    tic; % start clock
                    [stats, chains, diagnostics, info] = callbayes(engine, ...
                        'model'           , sprintf('%s_%s.txt', modelName, engine)   , ...
                        'data'            , data                                      , ...
                        'outputname'      , 'samples'                                 , ...
                        'init'            , generator                                 , ...
                        'datafilename'    , modelName                                 , ...
                        'initfilename'    , modelName                                 , ...
                        'scriptfilename'  , modelName                                 , ...
                        'logfilename'     , sprintf('tmp/%s', modelName)              , ...
                        'nchains'         , nChains                                   , ...
                        'nburnin'         , nBurnin                                   , ...
                        'nsamples'        , nSamples                                  , ...
                        'monitorparams'   , params                                    , ...
                        'thin'            , nThin                                     , ...
                        'workingdir'      , sprintf('tmp/%s', modelName)              , ...
                        'verbosity'       , 0                                         , ...
                        'saveoutput'      , true                                      , ...
                        'allowunderscores', true                                      , ...
                        'parallel'        , doParallel                                );
                    fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing
                    fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
                    if ~isfolder('storage')
                        !mkdir storage
                    end
                    save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
                    
                    % convergence of each parameter
                    disp('Convergence statistics:')
                    grtable(chains, 1.05)
                    
                    % basic descriptive statistics
                    disp('Descriptive statistics for all chains:')
                    codatable(chains);
                    
                end
                
                % convergence of each parameter
                disp('Convergence statistics:')
                grtable(chains, 1.05)
                
                % basic descriptive statistics
                disp('Descriptive statistics for all chains:')
                codatable(chains);
                
                % Savage-Dickey on delta
                
                % constants
                critical = 0;
                fontSize = 18;
                lineWidth = 2;
                markerSize = 10;
                posteriorColor = pantone.ClassicBlue;
                priorColor = pantone.Titanium;
                
                F = figure; clf; hold on;
                setFigure(F, [0.2 0.2 0.6 0.4], '');
                
                for idx = 1:2
                    
                    switch idx
                        case 1
                            binHi = 1; binLo = -1; binWidth = 0.005; binTick = 0.5;
                        case 2
                            binLo = -0.01; binHi = 0.01; binWidth = 0.001; binTick = 0.01;
                    end
                    binCenter = binLo:binWidth:binHi;
                    [~, criticalIdx] = min(abs(binCenter - critical));
                    
                    subplot(1, 2, idx); hold on;
                    
                    % axis
                    set(gca, ...
                        'ycolor'     , 'none'                , ...
                        'xlim'       , [binLo binHi]         , ...
                        'xtick'      , binLo:binTick:binHi   , ...
                        'box'        , 'off'                 , ...
                        'tickdir'    , 'out'                 , ...
                        'layer'      , 'top'                 , ...
                        'ticklength' , [0.02 0]              , ...
                        'layer'      , 'top'                 , ...
                        'fontsize'   , fontSize              );
                    set(gca, 'pos', get(gca, 'pos') + [-0.05 0.075 0 0]);
                    
                    % labels
                    xlabel('Difference in Rates', 'fontsize', fontSize+2);
                    
                    densityPrior = histcounts(chains.deltaPrior(:), ...
                        'binlimits'     ,[binLo-binWidth/2 binHi+binWidth/2]    , ...
                        'binwidth'      , binWidth         , ...
                        'normalization' , 'pdf'            );
                    plot(binCenter, densityPrior, '--', ...
                        'color'     , priorColor  , ...
                        'linewidth' , lineWidth   );
                    
                    densityPosterior = histcounts(chains.delta(:), ...
                        'binlimits'     , [binLo-binWidth/2 binHi+binWidth/2]    , ...
                        'binwidth'      , binWidth         , ...
                        'normalization' , 'pdf'            );
                    plot(binCenter, densityPosterior, '-', ...
                        'color'     , posteriorColor  , ...
                        'linewidth' , lineWidth   );
                    
                    
                    if idx == 1
                        legend('prior', 'posterior', ...
                            'box', 'off', ...
                            'fontsize', fontSize, ...
                            'location', 'northwest', ...
                            'autoupdate', 'off');
                    end
                    
                    BF = exp(log(densityPrior(criticalIdx)) - log(densityPosterior(criticalIdx)));
                    if idx == 2
                        if BF > 1
                            str = sprintf(' BF_{10} = %1.1f', BF);
                        else
                            str = sprintf(' BF_{01} = %1.1f', 1/BF);
                            
                        end
                        plot(binCenter(criticalIdx)*ones(1, 2), [densityPrior(criticalIdx) densityPosterior(criticalIdx)], 'k-');
                        text(binCenter(criticalIdx), mean([densityPrior(criticalIdx) densityPosterior(criticalIdx)]), str, ...
                            'fontsize', fontSize, ...
                            'horizontal', 'left');
                        
                        fprintf('%s: %s\n', d.roleNames{roleIdx}, str);
                    end
                    
                    plot(binCenter(criticalIdx), densityPrior(criticalIdx), 'o', ...
                        'linewidth'   , lineWidth, ...
                        'markeredgecolor',  priorColor, ...
                        'markerfacecolor', 'w', ...
                        'markersize', markerSize);
                    plot(binCenter(criticalIdx), densityPosterior(criticalIdx), 'o', ...
                        'linewidth'   , lineWidth, ...
                        'markeredgecolor',  posteriorColor, ...
                        'markerfacecolor', 'w', ...
                        'markersize', markerSize);
                    
                end
                % print
                if doPrint
                    warning off;
                    print(sprintf('figures/%s_%s_%s_SavageDickey.png', modelName, dataName, d.roleNames{roleIdx}), '-dpng');
                    print(sprintf('figures/%s_%s_%s_SavageDickey.pdf', modelName, dataName, d.roleNames{roleIdx}), '-dpdf');
                    print(sprintf('figures/%s_%s_%s_SavageDickey.eps', modelName, dataName, d.roleNames{roleIdx}), '-depsc');
                    warning on;
                end
                
            end
            
        case 'humansVsBots'
            
            winHtmp = []; propHnum = 0; propHden = 0;
            leagueList = unique(d.league);
            for leagueIdx = 1:length(leagueList)
                match = find(d.league == leagueList(leagueIdx));
                teamList = unique(d.team(match));
                humanTeamList = teamList(find(teamList >= 0));
                propHnum = propHnum + length(humanTeamList);
                propHden = propHden + length(teamList);
                for idx = 1:length(humanTeamList)
                    matchH = find(d.league == leagueList(leagueIdx) & d.team == humanTeamList(idx));
                    winHtmp = [winHtmp d.wins(matchH(1))./d.games(matchH(1))];
                end
            end
            winH = nanmean(winHtmp);
            propH = propHnum/propHden;
            % these two must average 0.5
            % propH*winH + (1-propH)*winB = 0.5 => winB = (0.5-propH*winH)/(1-propH)
            winB = (0.5-propH*winH)/(1-propH);
            fprintf('Human win %1.2f, bot win %1.2f\n', winH, winB);
            return; % gino, full story with no pics
            
        case 'valueCliffs'
            
            fontSize = 18;
            maxDepth = 60;
            dropPoint = 5;
            
            points = zeros(3e2, d.nRoles);
            denom = zeros(3e2, d.nRoles);
            % long lists seem to be provided from 2015 onwards
            yearList = 2015:2020;
            for yearIdx = length(yearList)
                T = readtable(sprintf('../data/FantasyPros_Fantasy_Football_Points_PPR_%d.csv', yearList(yearIdx)));
                T{:, 4} = strrep(T{:, 4}, 'DST', 'DS');
                for roleIdx = 1:d.nRoles
                    match = find(strcmp(T{:, 4}, d.roleNames{roleIdx}));
                    points(1:length(match), roleIdx) = points(1:length(match), roleIdx) + T{match, 5};
                    denom(1:length(match), roleIdx) = denom(1:length(match), roleIdx) + 1;
                end
            end
            avg = points./denom;
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.6 0.6], '');
            
            % composition axes
            set(gca, ...
                'xlim'          ,  [1 maxDepth]     , ...
                'xtick'         ,  [1 10:10:maxDepth]       , ...
                'ylim'          ,  [0 500]       , ...
                'ytick'         , 0:100:500       , ...
                'box'           ,  'off'                  , ...
                'tickdir'       ,  'out'                  , ...
                'ticklength'    ,  [0.01 0]               , ...
                'layer'         ,  'top'                  , ...
                'fontsize'      ,  fontSize                      );
            xlabel('Player Ranking Within Role', 'fontsize', fontSize+4);
            ylabel('Points for Season', 'fontsize', fontSize+4);
            set(gca, 'pos', get(gca, 'pos') + [0.01 0.05 0 0]);
            Raxes(gca, 0.03, 0.05);
            
            for roleIdx = 1:d.nRoles
                match = find(denom(:, roleIdx) == length(yearIdx));
                match = match(1:min(maxDepth, length(match)));
                val = avg(match, roleIdx);
                plot(1:length(match), val, 's-', ...
                    'color', roleColorsFull{roleIdx}, ...
                    'marker', roleMarkers{roleIdx}, ...
                    'markerfacecolor', roleColorsFull{roleIdx}, ...
                    'markeredgecolor', 'w', ...
                    'markersize', 10, ...
                    'linewidth', 2);
                text(length(match)+1, val(end), d.roleNames{roleIdx}, ...
                    'fontsize', fontSize-2);
                text(0, val(1), d.roleNames{roleIdx}, ...
                    'hor', 'right', ...
                    'fontsize', fontSize-2);
                fprintf('1st to %dth drop is %1.0f to %1.0f (%1.0f difference) for %s\n', dropPoint, val(1), val(dropPoint), val(dropPoint)-val(1), d.roleNames{roleIdx});
            end
            
        case 'roleOrder'
            
            maxSlots = 150;
            minPicks = 50;
            bounds = [25 75];
            winThreshold = [0.465 0.57];
            fontSize = 18;
            
            yearIdx = 2017;
            T = readtable(sprintf('../data/FantasyPros_Fantasy_Football_Points_PPR_%d.csv', yearIdx));
            T{:, 4} = strrep(T{:, 4}, 'DST', 'DS');
            
            m = zeros(d.nPlayers, maxSlots);
            value = nan(d.nPlayers, 1);
            pickMean = zeros(d.nPlayers, 1);
            for idx = 1:d.nPlayers
                match = find(d.player == (idx-1) & d.team >= 0 & d.auto ~= 1);
                if ~isempty(match)
                    grab = find(strcmp(d.playerNames{match(1)}, T{:, 2}));
                    if isempty(grab)
                        fprintf('Did not find %s\n', d.playerNames{match(1)});
                    elseif length(grab) > 1
                        fprintf('Found more than one %s\n', d.playerNames{match(1)});
                    else
                        value(idx) = T{grab, 5};
                    end
                    count = histcounts(d.pick(match), 0.5:maxSlots+1.5);
                    m(idx, 1:maxSlots) = count(1:end-1);
                    if sum(count) >= minPicks
                        pickMean(idx) = mean(d.pick(match));
                    else
                        pickMean(idx) = nan;
                    end
                    winProp{idx} = d.wins(match)./d.games(match);
                end
            end
            
            meanWins = cellfun(@mean, winProp);
            seWins = cellfun(@std, winProp)./sqrt(cellfun(@length, winProp));
            
            keep = find(~isnan(value) & ~isnan(pickMean));
            value = value(keep);
            pickMean = pickMean(keep);
            meanWins = meanWins(keep);
            seWins = seWins(keep);
            
            [~, sortIdx] = sort(pickMean, 'ascend');
            length(sortIdx)
            sortIdx = sortIdx(1:min(length(sortIdx), maxSlots));
            pickMean = pickMean(sortIdx);
            value = value(sortIdx);
            meanWins = meanWins(sortIdx);
            seWins = seWins(sortIdx);
            
            tmpRole = zeros(size(sortIdx));
            tmpPlayer = cell(size(sortIdx));
            for idx = 1:length(sortIdx)
                match = find(d.player == keep(sortIdx(idx))-1);
                tmpRole(idx) = d.role(match(1));
                tmpPlayer(idx) = d.playerNames(match(1));
            end
            
            % figure window
            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.55 0.6], '');
            
            for roleIdx = 1:d.nRoles
                subplot(2, 3, roleIdx); hold on;
                set(gca, ...
                    'xlim'          ,  [1 length(sortIdx)]         , ...
                    'xtick'         ,  [1 length(sortIdx)]   , ...
                    'ylim'          ,  [0 400]       , ...
                    'ytick'         ,  [0 100:100:400]         , ...
                    'box'           ,  'off'                  , ...
                    'tickdir'       ,  'out'                  , ...
                    'ticklength'    ,  [0.02 0]               , ...
                    'layer'         ,  'top'                  , ...
                    'fontsize'      ,  fontSize                      );
                set(gca, 'Position', get(gca, 'Position') + [0 0.05 0 0]);
                
                match = find(tmpRole == roleIdx);
                plot(pickMean(match), value(match), 'o', ...
                    'markersize', 8,  ...
                    'markerfacecolor', roleColorsFull{roleIdx}, ...
                    'color',  roleColorsFull{d.role(roleIdx)}, ...
                    'markeredgecolor', 'w');
                tmp = corrcoef(pickMean(match), value(match));
                text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), sprintf('r = %1.2f', tmp(1, 2)), ...
                    'fontsize', fontSize, ...
                    'vert', 'top', ...
                    'hor', 'right');
                
                text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), sprintf('%s', d.roleNames{roleIdx}), ...
                    'fontsize', fontSize, ...
                    'vert', 'top', ...
                    'hor', 'left');
                
                Raxes(gca, 0.01, 0.01);
                
            end
            
            [~, H(1)] = suplabel('Draft Position');
            [~, H(2)] = suplabel('Points Scored', 'y');
            set(H, 'fontsize', fontSize + 4);
            
    end % switch
    
    % Print figure
    if doPrint
        print(['figures/' aList{aListi} figAppend '.png'], '-dpng', '-r600');
        print([figDir aList{aListi} figAppend '.eps'], '-depsc');
    end
end