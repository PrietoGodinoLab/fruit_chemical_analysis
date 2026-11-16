%% rf_vs_pca_full.m
% Compare the top host-vs-non-host odours (random forest, 'mean' importance)
% with the top among-host odours (PCA, 'PC_score'), for BOTH feature-set
% criteria, in a single run:
%
%   CONSERVATIVE (stringent): RF trained only on odours measurable across all
%                datasets (excludes acids not seen in the SPME D. melanogaster set).
%   PERMISSIVE  : RF trained on all PDMS noni+pandan odours (admits the acids).
%
% For each criterion x species it:
%   1) reads the RF ranking and the PCA ranking,
%   2) keeps the odours common to both (matched by InChIKey),
%   3) takes the top N from each,
%   4) draws a WEIGHTED SANKEY (node height & ribbon width proportional to each
%      odour's share of its side's score), and
%   5) records the overlap for the significance figure.
% After each criterion it draws one OVERLAP-SIGNIFICANCE figure (noni + pandan),
% testing the overlap against a hypergeometric null.
%
% Outputs (vector PDF by default; set figFormat): 4 Sankeys + 2 significance figures.
% Base MATLAB only - no toolboxes. Put all input files (listed below) in the
% same folder as this script, or cd into that folder, then press Run.
%
% ---------------------------------------------------------------------------
% INPUT FILES REQUIRED IN THIS FOLDER
%   RF rankings (produced by rf_host_vs_nonhost.py):
%       noni_vs_others_RF_named_conservative.xlsx
%       noni_vs_others_RF_named_permissive.xlsx
%       pandan_vs_others_RF_named_conservative.xlsx
%       pandan_vs_others_RF_named_permissive.xlsx
%   PCA (among-host) rankings:
%       noni-pcscpre-inkey_named.xlsx
%       pandan-pcscore-inkey_named.xlsx
% ---------------------------------------------------------------------------

%% ============================ USER SETTINGS ==============================
N = 6;              % top-N per method (change to 7, 8, ... as desired)
saveFigs  = true;   % save every figure to a file next to this script
figFormat = 'pdf';  % 'pdf' (vector, editable in Illustrator) or 'svg'; use 'png' for a raster image

approaches = {'conservative','permissive'};   % labels used for files/titles

% Per species: the PCA file (shared by both criteria) and the RF file per criterion.
sp(1).label='Noni';   sp(1).pca='noni-pcscpre-inkey_named.xlsx';
sp(1).conservative='noni_vs_others_RF_named_conservative.xlsx';
sp(1).permissive ='noni_vs_others_RF_named_permissive.xlsx';

sp(2).label='Pandan'; sp(2).pca='pandan-pcscore-inkey_named.xlsx';
sp(2).conservative='pandan_vs_others_RF_named_conservative.xlsx';
sp(2).permissive ='pandan_vs_others_RF_named_permissive.xlsx';

% Column positions (1-based). RF files: InChIKey,name,mean ; PCA files: name,InChIKey,PC_score
rfInchiCol=1; rfNameCol=2; rfValCol=3;
pcaNameCol=1; pcaInchiCol=2; pcaValCol=3;

%% ============================ MAIN =======================================
for a = 1:numel(approaches)
    approach = approaches{a};
    clear stats
    for s = 1:numel(sp)
        rfFile  = sp(s).(approach);
        pcaFile = sp(s).pca;

        % --- read ---
        Trf=readAny(rfFile); Tpca=readAny(pcaFile);
        rfKey=normKey(colStr(Trf,rfInchiCol)); rfName=cleanName(colStr(Trf,rfNameCol)); rfVal=colNum(Trf,rfValCol);
        pcKey=normKey(colStr(Tpca,pcaInchiCol)); pcName=cleanName(colStr(Tpca,pcaNameCol)); pcVal=colNum(Tpca,pcaValCol);
        ok=rfKey~=""&~isnan(rfVal); rfKey=rfKey(ok); rfName=rfName(ok); rfVal=rfVal(ok);
        ok=pcKey~=""&~isnan(pcVal); pcKey=pcKey(ok); pcName=pcName(ok); pcVal=pcVal(ok);

        % --- keep common odours ---
        common=intersect(rfKey,pcKey);
        m=ismember(rfKey,common); rfKey=rfKey(m); rfName=rfName(m); rfVal=rfVal(m);
        m=ismember(pcKey,common); pcKey=pcKey(m); pcName=pcName(m); pcVal=pcVal(m);

        % --- top N of each ---
        n=min(N,min(numel(rfVal),numel(pcVal)));
        [~,o1]=sort(rfVal,'descend'); iR=o1(1:n);
        [~,o2]=sort(pcVal,'descend'); iP=o2(1:n);
        Lname=rfName(iR); Lkey=rfKey(iR); Lval=rfVal(iR);
        Rname=pcName(iP); Rkey=pcKey(iP); Rval=pcVal(iP);
        shared=intersect(Lkey,Rkey);

        fprintf('\n%s | %s: common=%d top%d shared=%d\n', sp(s).label, approach, numel(common), n, numel(shared));

        % --- weighted Sankey ---
        ttl=sprintf('%s (%s): RF vs PCA - top %d of %d common  (%d/%d shared)   |   width = share of score', ...
                    sp(s).label, approach, n, numel(common), numel(shared), n);
        figS=drawAlluvial(Lname,Lkey,Lval,Rname,Rkey,Rval,n,ttl,approach);
        if saveFigs, saveFig(figS, sprintf('%s_%s_sankey.%s', lower(sp(s).label), approach, figFormat)); end

        % --- stash significance stats ---
        Nc=numel(common); obs=numel(shared); E=n*n/Nc;
        pmf=arrayfun(@(x) hygePmf(Nc,n,n,x), 0:n);
        stats(s)=struct('label',sp(s).label,'N',Nc,'K',n,'obs',obs, ...
                        'E',E,'fold',obs/E,'pmf',pmf,'p',sum(pmf(obs+1:end)));
        fprintf('   E=%.2f fold=%.1fx p=%.4f\n', E, obs/E, stats(s).p);
    end

    % --- overlap-significance figure for this criterion (noni + pandan) ---
    figO=drawSignificance(stats, N, approach);
    if saveFigs, saveFig(figO, sprintf('%s_overlap_significance.%s', approach, figFormat)); end
end

%% ============================ LOCAL FUNCTIONS ============================
function fig = drawAlluvial(Lname,Lkey,Lval,Rname,Rkey,Rval,n,ttl,approach)
    x0=0; x1=1; bw=0.05; availH=9.0; gap=0.35;
    Lpos=stackPos(Lval,availH,gap); Rpos=stackPos(Rval,availH,gap);
    shared=intersect(Lkey,Rkey,'stable');
    cmap=lines(max(numel(shared),1));
    colOf=containers.Map('KeyType','char','ValueType','any');
    for k=1:numel(shared), colOf(char(shared(k)))=cmap(k,:); end
    grey=[0.82 0.82 0.82];
    fig=figure('Color','w','Position',[100 100 1000 640]);
    ax=axes(fig,'Position',[0.02 0.03 0.96 0.86]); hold(ax,'on');
    for i=1:n
        j=find(Rkey==Lkey(i),1);
        if ~isempty(j), ribbon(ax,x0+bw,Lpos(i,1),Lpos(i,2),x1-bw,Rpos(j,1),Rpos(j,2),colOf(char(Lkey(i)))); end
    end
    for i=1:n
        c=grey; if isKey(colOf,char(Lkey(i))), c=colOf(char(Lkey(i))); end
        drawBar(ax,x0,Lpos(i,:),bw,c);
        text(ax,x0-bw-0.03,mean(Lpos(i,:)),sprintf('%s  (%.3f)',Lname(i),Lval(i)), ...
            'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',8.3,'Interpreter','none');
    end
    for j=1:n
        c=grey; if isKey(colOf,char(Rkey(j))), c=colOf(char(Rkey(j))); end
        drawBar(ax,x1,Rpos(j,:),bw,c);
        text(ax,x1+bw+0.03,mean(Rpos(j,:)),sprintf('%s  (%.1f)',Rname(j),Rval(j)), ...
            'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8.3,'Interpreter','none');
    end
    text(ax,x0,availH+0.5,sprintf('Random Forest\n(mean, %s)',approach),'Interpreter','none', ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',11);
    text(ax,x1,availH+0.5,sprintf('PCA\n(PC_score)'),'Interpreter','none', ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',11);
    axis(ax,'off'); xlim(ax,[-1.35 2.25]); ylim(ax,[-0.4 availH+1.7]);
    title(ax,ttl,'FontSize',12,'FontWeight','bold'); hold(ax,'off');
end

function fig = drawSignificance(stats, N, approach)
    nD=numel(stats);
    fig=figure('Color','w','Position',[120 120 560*nD 470]);
    ymax=max(arrayfun(@(r) max(r.pmf), stats));
    for d=1:nD
        r=stats(d); ax=subplot(1,nD,d); hold(ax,'on'); xs=0:numel(r.pmf)-1;
        for x=xs
            if x>=r.obs, c=[0.75 0.22 0.17]; else, c=[0.72 0.77 0.82]; end
            rectangle(ax,'Position',[x-0.4,0,0.8,r.pmf(x+1)],'FaceColor',c,'EdgeColor','k','LineWidth',0.5);
        end
        plot(ax,[r.E r.E],[0 ymax*1.05],'--','Color',[0.17 0.24 0.31],'LineWidth',1.5);
        text(ax,r.E,ymax*0.99,sprintf(' chance E=%.2f',r.E),'Color',[0.17 0.24 0.31], ...
            'FontSize',9,'VerticalAlignment','top','Interpreter','none');
        text(ax,r.obs,max(r.pmf(r.obs+1)+ymax*0.10,ymax*0.5),sprintf('observed = %d',r.obs), ...
            'Color',[0.75 0.22 0.17],'FontWeight','bold','FontSize',10, ...
            'HorizontalAlignment','center','Interpreter','none');
        title(ax,sprintf('%s (%s):  %d/%d shared, pool N=%d\nfold %.1fx    p = %.3f', ...
            r.label,approach,r.obs,r.K,r.N,r.fold,r.p),'FontSize',12,'FontWeight','bold','Interpreter','none');
        xlabel(ax,'overlap in top list (odours in both)','Interpreter','none');
        if d==1, ylabel(ax,{'probability under the null','(random / independent rankings)'}); end
        xlim(ax,[-0.6 numel(r.pmf)-0.4]); ylim(ax,[0 ymax*1.08]); set(ax,'XTick',xs,'Box','off');
    end
    try, sgtitle(sprintf('Overlap significance - %s feature set  (red bars sum to p-value)',approach),'FontSize',12.5); catch, end
end

function pos = stackPos(vals,availH,gap)
    share=vals./sum(vals); usable=availH-(numel(vals)-1)*gap; h=share*usable;
    pos=zeros(numel(vals),2); y=availH;
    for i=1:numel(vals), pos(i,:)=[y,y-h(i)]; y=y-h(i)-gap; end
end
function drawBar(ax,xc,tb,bw,c)
    rectangle(ax,'Position',[xc-bw,tb(2),2*bw,tb(1)-tb(2)],'FaceColor',c,'EdgeColor','k','LineWidth',0.5);
end
function ribbon(ax,x0,yLt,yLb,x1,yRt,yRb,c)
    t=linspace(0,1,60); xs=x0+(x1-x0).*t; s=3*t.^2-2*t.^3;
    top=yLt+(yRt-yLt).*s; bot=yLb+(yRb-yLb).*s;
    fill(ax,[xs,fliplr(xs)],[top,fliplr(bot)],c,'FaceAlpha',0.55,'EdgeColor','none');
end
function p = hygePmf(N,K,n,x)
    if x<max(0,n-(N-K)) || x>min(n,K), p=0; return; end
    lnck=@(a,b) gammaln(a+1)-gammaln(b+1)-gammaln(a-b+1);
    p=exp(lnck(K,x)+lnck(N-K,n-x)-lnck(N,n));
end
function saveFig(fig,name)
    % Vector (pdf/eps/svg) stays fully editable in Illustrator, with transparency
    % preserved; png falls back to a 300-dpi raster.
    [~,~,ext]=fileparts(name); ext=lower(ext);
    isVec = any(strcmp(ext,{'.pdf','.eps','.svg'}));
    try
        if isVec
            exportgraphics(fig,name,'ContentType','vector','BackgroundColor','white');  % R2020a+ (svg: R2023b+)
        else
            exportgraphics(fig,name,'Resolution',300);
        end
    catch   % older MATLAB: use painters renderer for true vector output
        switch ext
            case '.pdf', print(fig,name,'-dpdf','-painters');
            case '.eps', print(fig,name,'-depsc','-painters');
            case '.svg', print(fig,name,'-dsvg','-painters');
            otherwise,   print(fig,name,'-dpng','-r300');
        end
    end
    fprintf('  saved %s\n',name);
end
function T = readAny(file)
    try, T=readtable(file,'VariableNamingRule','preserve'); catch, T=readtable(file); end
end
function s = colStr(T,idx)
    c=T{:,idx};
    if isstring(c), s=c; elseif iscell(c), s=string(c); else, s=string(c); end
    s=s(:);
end
function v = colNum(T,idx)
    c=T{:,idx};
    if isnumeric(c), v=double(c); else, v=str2double(string(c)); end
    v=v(:);
end
function k = normKey(s), k=upper(strtrim(strrep(s,char(160),' '))); k(ismissing(k))=""; end
function nm = cleanName(s), nm=strtrim(strrep(s,char(160),' ')); nm(ismissing(nm))=""; end
