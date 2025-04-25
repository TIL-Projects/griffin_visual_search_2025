et = read.csv('aois_analyze_all_out_working.csv',as.is=T) 
names(et) = tolower(names(et))

b = read.csv('b.csv',as.is=T)
names(b) = tolower(names(b))

# get the eeg variables
m = read.csv('bioCorePull_ET_ages.csv',as.is=T)
m = m[,c('visit_individual','visit_time_period',grep('^(eeg|erp)',names(m),value=T))]
m = m[m$visit_time_period=='t1',]

manifest = read.csv('manifest.csv',as.is=T)
names(manifest) = tolower(names(manifest))

# combine all the files together
x = et
dim(x)
files_in_manifest = sort(unique(manifest$edf))
files_in_et = sort(unique(et$filebase))
files_in_manifest_but_not_et = setdiff(files_in_manifest,files_in_et)

x = merge(x,manifest,by.x='filebase',by.y='edf',all.x=T,suffixes=c('','.manifest'))
table(x$visit,useNA='ifany')
dim(x)
x = merge(x,b[,setdiff(names(b),grep('^et',names(b),value=T))],by.x='subject',by.y='visit_individual',all.x=T,suffixes=c('','.b'))
dim(x)
x = merge(x,m,by.x='subject',by.y='visit_individual',all.x=T,suffixe=c('','.m'))
dim(x)

table(visit=x$visit,visit_time_period=x$visit_time_period,useNA='ifany')
dim(x)

# remove t2 and t3 data
x = x[!is.na(x$visit),]

# make variables consistent
x$static_nn_mean = as.numeric(x$static_nn_mean)

x0 = x

# remove invalid data
x = x0
roi_dur_fields = grep('^valid_dur',names(x),value=T)
roi_dur_fields_valid = setdiff(roi_dur_fields,c('valid_dur_offscreen','valid_dur_undefined'))
x$valid_dur_all = apply(x[,roi_dur_fields_valid],1,sum)
x$pvalid = x$valid_dur_all/x$t1_t0

x$valid_time = x$pvalid>=.5
x$valid_calerr = x$calmethod!='invalid' & x$calerr<=100
x$valid_center_looking = x$latency_roi_1=='visualsearch_center'
x$valid = x$valid_time & x$valid_calerr & x$valid_center_looking

# transform visit information into long form
library(tidyverse)
x$key = paste(x$filebase,x$sourcefile)
table(duplicated(x$key))
latency_fields = grep('^latency',names(x),value=T)
yb = x[,c('key','dx','valid',latency_fields)]
yl = pivot_longer(yb,latency_fields,
    names_to=c('visit_field','visit_num'),
    names_pattern='^latency_([a-z0-9]+)_(\\d+)$',
    values_to='value',
    values_transform=list(value=as.character),
    values_drop_na=T
    )
yw = pivot_wider(yl,
    id_cols=c(key,dx,valid,visit_num),
    names_from='visit_field',
    names_prefix='visit_',
    values_from='value'
    )
LOST_DATA_INTERPOLATE_MAX_MS = 80
MAX_TIME_MS = 20000
SAMPLE_TIME_MS = 1/500*1000
numeric_fields = c('visit_num','visit_t0','visit_t1')
yw[,numeric_fields] = sapply(yw[,numeric_fields],as.numeric)

# aggregate data and strip out "non region" areas for analyses
yx = yw
yx[yx$visit_roi %in% c('visualsearch_eyes','visualsearch_head','visualsearch_mouth'),'visit_roi']='visualsearch_face'
yx = yx[!(yx$visit_roi %in% c('','visualsearch_background','offscreen')),]

# strip out center looks beyond the first visit
yx = yx[yx$visit_roi!='visualsearch_center' | yx$visit_num==1,]

# make sure the t0 is 0 at min
yx[yx$visit_t0<0,'visit_t0'] = 0


y = yx %>% mutate(
    visit_roi = gsub('faceoutline','outline',gsub('visualsearch_','',visit_roi)),
    visit_t1_t0 = visit_t1-visit_t0+SAMPLE_TIME_MS, # sample time added as correction 
    next_dt=if_else( # time to get to the next object
        lead(key)==key,
        pmax(0,lead(visit_t0)-visit_t1-SAMPLE_TIME_MS),
        pmax(0,as.double(MAX_TIME_MS-visit_t1-SAMPLE_TIME_MS/2))
        ), 
    first_roi_in_train=lag(key)!=key | lag(visit_roi)!=visit_roi,
    NULL)
y[1,'first_roi_in_train'] = T

# compress all repeated visit trains into a single continguous visit and get summary stats
y2 = y %>% 
    group_by(key) %>% 
    mutate(
        visit_num_org=visit_num,
        visit_num=cumsum(first_roi_in_train),
        next_dt0=next_dt,
        NULL
        ) %>%
    mutate(
        trial_max_visit_num=max(visit_num),
        NULL 
        ) %>%
    group_by(key,visit_num) %>% 
    summarise( .groups='keep',
        dx=first(dx),
        valid=first(valid),
        visit_num=first(visit_num),
        trial_max_visit_num=first(trial_max_visit_num),
        visit_num_org0=first(visit_num_org),
        visit_num_org1=last(visit_num_org),
        visit_roi=first(visit_roi),
        visit_t0=min(visit_t0),
        visit_t1=max(visit_t1),
        visit_t1_t0=visit_t1-visit_t0+SAMPLE_TIME_MS, # adjusted for sample uncertainty
        next_dt=last(next_dt0),
        NULL
        )

# compute prior transition times leading up to the current visit
y3 = y2 %>% 
    group_by(key) %>% 
    mutate(
        prior_cdt=cumsum(next_dt)-next_dt, # prior cumulative transition time
        prior_mdt=ifelse(visit_num>1,prior_cdt/(visit_num-1),NA), # avg transition time per transition
        #prior_cnt=visit_num-1, # prior number of transition // computed just for convenience
        prior_cobjtime=cumsum(visit_t1_t0)-visit_t1_t0, # prior time spent in objects (including gaps)
        prior_mobjtime=ifelse(visit_num>1,prior_cobjtime/(visit_num-1),NA), # avg time per obj
        #prior_cnobjs=visit_num-1, # prior number of objects = prior # transitions // convenience
        NULL
        )
rois = unique(y3$visit_roi)
for (roi in rois) {
    tvar = paste('ctime_',roi,sep='')
    nvar = paste('cn_',roi,sep='')
    y3 = y3 %>% 
    mutate(
        !!tvar := cumsum(visit_t1_t0*(visit_roi==roi)),
        !!nvar := cumsum(visit_roi==roi),
    )
}

# aggregate statistics (at top level)
y_all = y3 %>% group_by(key) %>% 
    summarise(
        dx=first(dx),
        valid=first(valid),
        num_visits = max(visit_num),
        last_visit_num0 = max(visit_num_org1),
        visited_center = any(visit_roi=='center'),
        num_regions_visited = length(unique(visit_roi)),
        NULL)

# info about overall region based behavior
y_region = y3 %>% 
    group_by(key,visit_roi) %>% 
    summarise(.groups='keep',
        dx=first(dx),
        valid=first(valid),
        num_looks=n(),
        mean_look_time=mean(visit_t1_t0,na.rm=T),
        med_look_time=median(visit_t1_t0,na.rm=T),
        sd_look_time=sd(visit_t1_t0,na.rm=T),
        first_look_visit_num=first(visit_num),
        first_look_t0=first(visit_t0),
        first_look_t1=first(visit_t1),
        first_look_t1_t0=first(visit_t1_t0),
        last_look_visit_num=last(visit_num),
        last_look_t0=last(visit_t0),
        last_look_t1=last(visit_t1),
        last_look_t1_t0=last(visit_t1_t0),
        look_expanse_visit_num=last(visit_num)-first(visit_num)+SAMPLE_TIME_MS,
        look_expanse_t0_t0=last(visit_t0)-first(visit_t0)+SAMPLE_TIME_MS,
        look_expanse_t0_t1=last(visit_t1)-first(visit_t0)+SAMPLE_TIME_MS,
        NULL)
y_regionw = y_region %>% pivot_wider(
   id_cols=c(key,dx,valid),
   names_from=visit_roi,
   values_from=num_looks:look_expanse_t0_t1,
   names_sep='_',
   values_fill=NA)


# info about timing of unique visits/first visits
y_uniq = y3 %>% filter(visit_roi!='center') %>% distinct(key,visit_roi,.keep_all=T) %>%
    group_by(key) %>% 
    mutate(
        unique_visit_num=row_number(),
        next_unique_dt=lead(visit_t0)-visit_t0,
        NULL)

# invert/transpose y_uniq so we can look at latencies of orienting to each type of object, NAs if not looked at
y_roi = y_uniq %>%
    pivot_wider(
        id_cols=c(key),
        names_from=visit_roi,
        values_from=c(
            unique_visit_num,visit_num,
            visit_t0:last_col()),
        names_sep='_',
        values_fill=NA)
y_uvisit = y_uniq %>%
    pivot_wider(
        id_cols=c(key),
        names_from=unique_visit_num,
        names_prefix='uv',
        values_from=c(
            visit_roi,
            visit_num,
            visit_t0:last_col()),
        names_sep='_',
        values_fill=NA) 
for (roi in rois) {
    for (uv in 1:5) {
        is_roi_name = paste('uv',uv,'_is_',roi,sep='')
        saw_roi_name = paste('uv',uv,'_saw_',roi,sep='')
        last_saw_roi_name = paste('uv',uv-1,'_saw_',roi,sep='')
        visit_roi_var = paste('visit_roi_uv',uv,sep='')
        y_uvisit = y_uvisit %>% 
            mutate(
                !!is_roi_name := !!sym(visit_roi_var)==roi,
                !!saw_roi_name := ifelse(
                    uv>1, 
                    (!!sym(is_roi_name)) | (!!sym(last_saw_roi_name)), 
                    !!sym(is_roi_name))
            ) 
    }
}

library(psych)
# create aggrgeate varaibles for analyses
my_mean = function(x) {
    mean(x[is.finite(x)])
}
agg_roi = y_uniq %>%
    filter(valid) %>%  
    mutate(
        filebase=gsub(' .*$','',key),
        sourcefile=gsub('^[^ ]+ ','',key),
        NULL
        ) %>% 
    group_by(filebase,visit_roi) %>% 
    summarise(.groups='keep',
        key=first(key),
        valid=first(valid),
        dx=first(dx),
        ntrials=n(),
        visit_t0=mean(visit_t0),
        prior_cobjtime=mean(prior_cobjtime),
        prior_mobjtime=my_mean(prior_cobjtime/(visit_num-1)),
        unique_visit_num=mean(unique_visit_num),
        visit_num=mean(visit_num))
agg_w = agg_roi %>% pivot_wider(
    id_cols=c(filebase),
    id_expand=F,
    names_prefix='agg',
    names_from=visit_roi,
    values_from=c(ntrials,visit_t0:last_col()),
    names_sep='_',
    values_fill=NA
    )

dim(x)
xm = merge(x,y_roi,by='key',all.x=T,suffixes=c('','.yroi'))
dim(xm)
xm = merge(xm,y_uvisit,by='key',all.x=T,suffixes=c('','.yuvisit'))
dim(xm)
xm = merge(xm,y_regionw,by='key',all.x=T,suffixes=c('','.yregionw'))
dim(xm)
xm = merge(xm,y_all,by='key',all.x=T,suffixes=c('','.yall'))
dim(xm)

#xm = merge(xm,agg_w[,c('filebase',setdiff(names(agg_w),names(xm)))],by='filebase',all.x=T)
#dim(xm)


x = xm

# compute variables of interest
x$fix_cnt_visualsearch_oem = 
    x$fix_cnt_visualsearch_head+ 
    x$fix_cnt_visualsearch_eyes+ 
    x$fix_cnt_visualsearch_mouth 
x$fix_dur_visualsearch_oem = 
    x$fix_dur_visualsearch_head+ 
    x$fix_dur_visualsearch_eyes+ 
    x$fix_dur_visualsearch_mouth 

x$looked_at_oem = 0+(x$fix_cnt_visualsearch_oem>=1)

oem_rlat_fields = c(
    'rlat_t0_1_visualsearch_head',
    'rlat_t0_1_visualsearch_eyes',
    'rlat_t0_1_visualsearch_mouth')

minfinite = function(x) {
    r = min(x,na.rm=T)
    if (!is.finite(r)) {
        return(NA)
    } else {
        return(r)
    }
}

sumfinite = function(x) {
    r = sum(x,na.rm=T)
    if (!is.finite(r)) {
        return(NA)
    } else {
        return(r)
    }
}

x$time_to_oem = apply(x[,oem_rlat_fields],1,minfinite)
x$time_to_em = apply(x[,setdiff(oem_rlat_fields,'rlat_t0_1_visualsearch_head')],1,minfinite)
x$time_to_outer = x$rlat_t0_1_visualsearch_head
x$time_to_eyes = x$rlat_t0_1_visualsearch_eyes
x$time_to_mouth = x$rlat_t0_1_visualsearch_mouth
x$tem = 
    x$valid_dur_visualsearch_eyes + 
    x$valid_dur_visualsearch_mouth  
x$toem = 
    x$valid_dur_visualsearch_head + 
    x$tem
x$poem = x$toem/x$valid_dur_all
x$pem = x$tem/x$valid_dur_all

# exploration - number of images explored
non_social_obj = c('bird','car','tech','faceoutline')
fix_cnt_fields_non_social = sprintf('fix_cnt_visualsearch_%s',non_social_obj)
looked_at_fields_non_social = gsub('fix_cnt_visualsearch','looked_at',fix_cnt_fields_non_social)
x[,looked_at_fields_non_social] = 0+(x[,fix_cnt_fields_non_social]>=1)
x$num_nonsocial_looked_at = apply(x[,looked_at_fields_non_social],1,sumfinite)
x$num_social_looked_at = 0+(x[,'fix_cnt_visualsearch_oem']>=1)

# perseveration - total fixation duration per image explored
fix_dur_fields_non_social = sprintf('fix_dur_visualsearch_%s',non_social_obj)
x$fix_dur_non_social = apply(x[,fix_dur_fields_non_social],1,sumfinite)
x$perseveration_non_social = x$fix_dur_non_social / x$num_nonsocial_looked_at
x$perseveration_social = x$fix_dur_visualsearch_oem / x$num_social_looked_at

# detail orientation - num_fixations / objects explored
x$fix_cnt_non_social = apply(x[,fix_cnt_fields_non_social],1,sumfinite)
x$detail_oriented_non_social = x$fix_cnt_non_social / x$num_nonsocial_looked_at
x$detail_oriented_social = x[,'fix_cnt_visualsearch_oem'] / x$num_social_looked_at

# convert bad values to NAs
fields_to_na = c('perseveration_non_social','perseveration_social','detail_oriented_non_social','detail_oriented_social','time_to_oem','time_to_eyes','time_to_mouth')
for (field in fields_to_na) {
    x[!is.finite(x[,field]),field] = NA
}

# normalize variables by pvalid
fields_to_normalize = c('num_nonsocial_looked_at','num_social_looked_at','perseveration_non_social','perseveration_social','detail_oriented_non_social','detail_oriented_social')
for (field in fields_to_normalize) {
    field_new = paste(field,'_norm',sep='')
    x[,field_new] = x[,field]/x$pvalid
}

x = x[x$valid,]
x1 = x

# aggregate data
x = x1
x$count = 1
super_aggregate = function (d,by) {
    var_types = unlist(lapply(d,class))
    var_names = names(var_types)

    idx = var_types=='Date'
    d[,idx] = lapply(d[,idx],as.character) # for simplicity turn dates into characters
    idx = var_types=='logical' | var_types == 'integer'
    d[,idx] = 0+d[,idx] # for simplicity turn logicals and integers iunto numeric
    var_types = unlist(lapply(d,class))
    var_names = names(var_types)

    idx = var_types=='character'
    a = aggregate(d[,idx],by,function(x) {paste(unique((x)),collapse='|')})

    idx = var_types=='numeric'
    b = aggregate(d[,idx],by,function(x) {mean(x,na.rm=T)})
    
    counts = aggregate(d[,1],by,function(x) {length(x)})
    names(counts)[length(by)+1] = 'count'

    y = cbind(counts,a[,-(1:length(by))],b[,-(1:length(by))])
    y = y[,c(setdiff(names(y),names(d)),names(d))]
    return(y);
}

bylist = list(filebase=x$subject,timepoint=x$visit)
y = super_aggregate(x,bylist)
x2 = y

# remove aggregated data that fails trial inclusion criteria
x = x2
min_num_trials = .25 * 12
x$valid_trials = x$count>=min_num_trials
table(x$valid_trials,x$dx)
x = x[x$valid_trials,] # only kids who have enough trials to analyze

library(psych)
x3 = x

#-------------- start statistical analyses
x = x3
library(car)
library(emmeans)
options(contrasts = c("contr.helmert", "contr.poly"));
x$dx = as.factor(x$dx)
contrasts(x$dx) = contr.sum
x$sex = 'unk'
x[x$is_male==1,'sex'] = 'male'
x[x$is_male==0,'sex'] = 'female'
x$sex = as.factor(x$sex)
contrasts(x$sex) = contr.sum

y = x
describeBy(y[,c('age_years','is_male','visit_ados2_sa_severity','visit_ados2_rrb_severity','visit_ados2_comparison_score_ss','imported_nepsy_calc_mf_scaled_score','imported_srs_t_sci','imported_srs_t_score.t1','diagnosis_summary_verbal_best_iq','diagnosis_summary_nonverbal_best_iq','diagnosis_summary_full_scale_best_iq')],group=list(dx=y$dx))

library(lme4)
adj='none'
x = x[x$timepoint=='t1',]
x$subject = as.factor(as.character(x$subject))
stats_of = function(varname) {
    x$var = x[,varname]
    hist(x$var)

    print('--------------------------')
    print(varname)
    print('--------------------------')

    print(describeBy(x$var,group=list(dx=x$dx)))

    print('***************************') 
    print('**** without controls')
    #m = lmer(var~dx+(1|subject),data=x)
    m = glm(log(var)~dx,data=x)
    a=Anova(m,3,test='F'); s=summary(m)
    print(a)
    print(s)
    m.s=emmeans(m,~dx,type='response');p=pairs(m.s,adj=adj)
    #m.s=emmeans(m,~sex,type='response');p=pairs(m.s)
    print(m.s)
    print(p)
    p.s=summary(p)
    #p.s$estimate = p.s$estimate/sigmaHat(m)
    #print(p.s)

    print('***************************') 
    print('**** with controls')
    m = glm(log(var)~dx+t1iq+age_years+pvalid,data=x)
    a=Anova(m,3,test='F'); s=summary(m)
    print(a)
    print(s)
    m.s=emmeans(m,~dx,type='response');p=pairs(m.s,adj=adj)
    print(m.s)
    print(p)
    p.s=summary(p)
    #p.s$estimate = p.s$estimate/sigmaHat(m)
    #print(p.s)

    print('***************************') 
    print('**** with controls and sex')
    m = glm(log(var)~dx*sex+t1iq+age_years+pvalid,data=x)
    a=Anova(m,3,test='F'); s=summary(m)
    print(a)
    print(s)
    m.s=emmeans(m,~dx,type='response');p=pairs(m.s,adj=adj)
    print(m.s)
    print(p)
    p.s=summary(p)
    m.s=emmeans(m,~sex,type='response');p=pairs(m.s,adj=adj)
    print(m.s)
    print(p)
    p.s=summary(p)
    m.s=emmeans(m,~dx*sex,type='response');p=pairs(m.s,adj=adj)
    print(m.s)
    print(p)
    p.s=summary(p)
}

stats_of('visit_t0_face')
stats_of('visit_num_face')
stats_of('prior_mobjtime_face')
