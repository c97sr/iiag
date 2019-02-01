# Formats output from the functions before printing it out to the file
# Since the three methods return different number of variables
# and the order of the columns in each variable is not the same
# even when the columns exist, we format all result set variables
# to common format
# Credit: Sasikiran Kandula, Columbia University Mailman School of Public Health,
#         Department of Environmental Health Sciences

buildOutput<-function(res_sub, format='long'){
  required_order=NULL
  if(format == 'long'){
    required_order = c("fc_start","time","week","Est","Est_sd","S","S_sd","I","I_sd","L","L_sd","D","D_sd","R0max","R0max_sd","R0min","R0min_sd");
  }else{
    required_order = c("fc_start","time","week","Est","Est_sd");
  }
  
  out=subset(res_sub, select=required_order[1]);#init with fc_start which all result sets have
  
  for(i in 2:length(required_order)){
    if(required_order[i] %in% colnames(res_sub)){
      #index into res using name of the column and bind it to existing output, round to 0/2 digits
      if(required_order[i] %in% c("S","I","Est","E", "S_sd","I_sd","Est_sd","E_sd")){
        d=0
      }else{
        d=2
      }
        
      out = cbind(out, round(res_sub[,required_order[i]], digits=d))
      
    } else{# a column of -1s
      out = cbind(out, matrix(-1, nrow(res_sub),1))
      
    }
  }
  
  colnames(out)=required_order
  out
}

# Converts time column in the result set from 'number of days since start of year' to 'date'
# Calculate dates by adding number of days (the current value in time column) to 01-01 of current year
# The current year is calculated by extracting the year part of start date which is the season start
# Output has an additional column next to 'time'
convertDaysToDate<-function(res, start_date){
  #get *.time column which has days i.e. clim_start.
  #convert this into date using week_start
  time_index = grep("*time", colnames(res));
  if(length(time_index) != 1){#Could not find the index, return result as is
    out = res;
  }
  
  else{
    time_index = as.numeric(time_index[1]);
    time_column = NULL;
    
    for(id in 1:nrow(res))
      time_column = append(time_column, as.Date(res[id, time_index]-1, origin=paste(1900+as.POSIXlt(start_date)$year,"01","01",sep="-")));

    time_column = data.frame(time_column);
    
    out=cbind(res[,1:time_index], time_column, res[,(time_index+1):ncol(res)])
    
  }
  
  out
}

# This function takes as argument 
# a matrix to be written out, 
# a target (filename if csv, tablename if db), 
# and a type (metrics or output)
writeOutput <- function(output, targetName, type){
  if(output_conn_type == "csv"){#Write to metrics file
    write.table(output, targetName, append=T, col.names=F, row.names=F, sep=",")
  }
  
  else if(output_conn_type == "db"){#Write to DB
    
    if(type == "metrics"){
      field_types = metrics_field_types
    }
    
    else if(type == "output"){
      field_types = output_field_types
    }
    
    
    if(dbExistsTable(con, targetName))
      dbWriteTable(con, targetName, as.data.frame(output), append = T, row.names = F)
    else
      dbWriteTable(con, targetName, as.data.frame(output), field.types = field_types, row.names = F)
  }
}

# Takes an array of newI and a baseline
# Onset = index into values where value > baseline for at least 3 consecutive weeks
# End = index into values where the period of elevate flu incidence ends i.e
# all indices between onset and end(both inclusive) satisy baseline condition
# Duration = end - offset + 1
# All NA if no 3 week window was found
# Need to add week offset to return values to get the actual weeks
# Also assumes no double-peaks i.e. at most one window satisfies condition
findOnset <- function(values, baseline){
  onset = end = duration = NA
  
  above = which(values > baseline);
  
  if(length(above) > 2){
    for(i in 1:length(above)){
      if((above[i]+1) %in% above && (above[i]+2) %in% above){
        onset = above[i] + wk_start - 1;
        break;
      }
    }
    
    for(i in length(above):1){
      if((above[i]-1) %in% above && (above[i]-2) %in% above){
        end = above[i] + wk_start - 1;
        break;
      }
    }
    
    if(!is.na(onset) && !is.na(end)){
      duration = end - onset + 1
    }
    
    else{
      onset = end = duration = NA
    }
  }
  
  out = list(onset =onset, end = end, duration = duration)
  #print(out)
}

formatDist = function(onsetsDist, peakWeeksDist, peakIntensitiesDist, nextILIDist){
  onset = cbind('onset', onsetsDist)
  pw = cbind('pw', peakWeeksDist)
  week = rbind(onset, pw)
  
  pi = cbind('pi', peakIntensitiesDist)
  for(i in 2:ncol(nextILIDist)){
    type = paste('next_week_', (i-1), sep="")
    pi = rbind(pi, cbind(type, nextILIDist[, 1], nextILIDist[, i]))
  }
  
  return(rbind(week, pi))
}

formatDistNEW = function(onsets3Dist, onsets4Dist, onsets5Dist, onsets6Dist,
                         peakWeeksDist, peakIntensitiesDist, nextILIDist){
  onset3 = cbind('onset3', onsets3Dist)
  onset4 = cbind('onset4', onsets4Dist)
  onset5 = cbind('onset5', onsets5Dist)
  onset6 = cbind('onset6', onsets6Dist)
  pw = cbind('pw', peakWeeksDist)
  week = rbind(onset3, onset4, onset5, onset6, pw)
  
  pi = cbind('pi', peakIntensitiesDist)
  for(i in 2:ncol(nextILIDist)){
    type = paste('next_week_', (i-1), sep="")
    pi = rbind(pi, cbind(type, nextILIDist[, 1], nextILIDist[, i]))
  }
  
  return(rbind(week, pi))
}

formatDistENDS = function(onsetsDist, endsDist, dursDist, peakWeeksDist, peakIntensitiesDist, nextILIDist){
  onset = cbind('onset', onsetsDist)
  end = cbind('end', endsDist)
  pw = cbind('pw', peakWeeksDist)
  week = rbind(onset, end, pw)
  
  dur = cbind('dur', dursDist)
  pi = cbind('pi', peakIntensitiesDist)
  for(i in 2:ncol(nextILIDist)){
    type = paste('next_week_', (i-1), sep="")
    pi = rbind(pi, cbind(type, nextILIDist[, 1], nextILIDist[, i]))
  }
  
  return(rbind(week, dur, pi))
}