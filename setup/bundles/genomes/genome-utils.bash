function get_job {

  tool_id=$1

  jobs="$(chipster -o json job list)"
        
  for j in $(seq $(echo "$jobs" | jq length)); do     
    job_tool_id="$(echo $jobs | jq .[$j-1].toolId -r)"
    job_id="$(echo $jobs | jq .[$j-1].jobId -r)"
    
    if [ "$tool_id" = "$job_tool_id" ]; then
      break
    else
      job_id=""
    fi
  done
  
  if [ -z "$job_id" ]; then
    >&2 echo "job $tool_id not found"
    exit 1
  fi
  
  echo $job_id
}

# Define sensible packages

# An own package will be create for each index. This job creates only the extracted folder which can be used directly in development machines. 
# The packaging itself will happen in another job. But we will create here 
# a file which defines the appropriate packages because here we have the necessary information of genomes and versions available. 
function add_file_to_package {
  dir="$1"
  file="$2"
  package="$3"
  
  long_package=$(echo $dir | tr "/" "_")
      
  if [ -n "$package" ]; then
    long_package="${long_package}_${package}"
  fi
        
  echo "$long_package	$dir/$file" >> $packages_file
}

function download_result_datasets {
  tool_id="$1"
  dir="$2"
  package="$3"
  
  if [ -n "$dir" ]; then
    mkdir -p $dir  
    pushd $dir
  fi
  
  job_id=$(get_job $tool_id)
  
  datasets="$(chipster -o json dataset list)"
        
  for j in $(seq $(echo "$datasets" | jq length)); do

    source_job="$(echo $datasets | jq .[$j-1].sourceJob -r)"
    dataset_id="$(echo $datasets | jq .[$j-1].datasetId -r)"
    dataset_name="$(echo $datasets | jq .[$j-1].name -r)"
    
    if [ "$source_job" = "$job_id" ]; then
      echo "download $dataset_name"
      
      chipster --quiet dataset download --file "$dataset_name" $dataset_id
      
      if [ "$package" != "none" ]; then
        add_file_to_package "$dir" "$dataset_name" "$package"
      fi
    fi
  done
  
  if [ -n "$dir" ]; then
    popd
  fi
}

# Create symlink and define a package for it
function create_symlink {
  target="$1"
  link_dir="$2"
  link_name="$3"
  package="$4"
  
  ln -s $target $link_dir/$link_name
  
  add_file_to_package "$link_dir" "$link_name" "$package"
}

function wait_for_jobs {

  genomes_path="$1"

  echo "wait jobs to finish at $(date --iso-8601=seconds)"

  ready=0

  while [ $ready = 0 ]; do
    ready=1
    
    echo ""
    echo "$(date --iso-8601=seconds)"
    echo "WAIT. 	SCHED. 	RUN. 	SPECIES 	VERSION 	RELEASE"
    echo "---------------------------------------------------------------------------"
    while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT; do

        session_name="$JOB_NAME-$BUILD_NUMBER-$SPECIES-$VERSION-$RELEASE"
        
        # don't stop even if some status checks fail
        set +e
      
        chipster --quiet session open $session_name
      
        jobs="$(chipster --quiet job list)"
    
        new=$(echo "$jobs" | grep NEW | wc -l)
        waiting=$(echo "$jobs" | grep WAITING | wc -l)
        scheduled=$(echo "$jobs" | grep SCHEDULED | wc -l)
        running=$(echo "$jobs" | grep RUNNING | wc -l)
    
        if [ $new != 0 ] || [ $waiting != 0 ] || [ $scheduled != 0 ] || [ $running != 0 ]; then
          ready=0
        
          echo "$waiting 	$scheduled 	$running 	$SPECIES 	$VERSION 	$RELEASE"
        fi
        
        set -e
        
    done < $genomes_path
    
    if [ $ready = 0 ]; then
      sleep 30
    fi
  done
}
