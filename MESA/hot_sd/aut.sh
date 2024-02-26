
#!/bin/bash

function modify_inlist_to_TAMS {
    echo -e "\nZbase=$new_Zbase,M=$new_mass"

    sed -i "s#save_model_filename = .*#save_model_filename = './mod_TAMS/TAMS_$1-$2.mod'#" inlist_to_TAMS

    sed -i "s/initial_mass = .*/initial_mass = $1/" inlist_to_TAMS
    sed -i "s/Zbase = .*/Zbase = $2/" inlist_to_TAMS
    sed -i "s/initial_z = .*/initial_z = $2/" inlist_to_TAMS

    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
    do_one inlist_to_TAMS_header
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
}

function modify_inlist_to_ZACHeB {
    echo -e "\nZbase=$new_Zbase,M=$new_mass"

    sed -i "s#load_model_filename = .*#load_model_filename = './mod_TAMS/TAMS_$1-$2.mod'#" inlist_to_ZACHeB
    sed -i "s#save_model_filename = .*#save_model_filename = './mod_pre_zahb/pre_ZAHB_$1-$2.mod'#" inlist_to_ZACHeB

    sed -i "s/initial_mass = .*/initial_mass = $1/" inlist_to_ZACHeB
    sed -i "s/Zbase = .*/Zbase = $2/" inlist_to_ZACHeB
    sed -i "s/initial_z = .*/initial_z = $2/" inlist_to_ZACHeB
    
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
    do_one inlist_to_ZACHeB_header
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

    cp LOGS_POSTMS/history.data data_pre/$1-$2.data
    rm -f ./LOGS_POSTMS/*
}

function modify_inlist_remove_env {
    echo -e "\nZbase=$new_Zbase,M=$new_mass,env_mass=$new_env_mass"

    sed -i "s#load_model_filename = .*#load_model_filename = './mod_pre_zahb/pre_ZAHB_$1-$2.mod'#" inlist_remove_env
    sed -i "s#save_model_filename = .*#save_model_filename = './mod_env/ZAHB_$1-$3-$2.mod'#" inlist_remove_env

    sed -i "s/initial_mass = .*/initial_mass = $1/" inlist_remove_env
    sed -i "s/Zbase = .*/Zbase = $2/" inlist_remove_env
    sed -i "s/initial_z = .*/initial_z = $2/" inlist_remove_env
    sed -i "s/extra_mass_retained_by_remove_H_env = .*/extra_mass_retained_by_remove_H_env = $3/" inlist_remove_env

    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
    do_one inlist_remove_env_header
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

    cp LOGS_env/history.data data_remove_env/$1-$3-$2.data
    rm -f ./LOGS_env/*
}

function modify_inlist_hb {
    echo -e "\nZbase=$new_Zbase,M=$new_mass,env_mass=$new_env_mass,Y_surf=$Y_surf"
    
    result=$(python3 aut_q.py $new_mass $new_env_mass $new_Zbase)
    echo "Result from Python script: $result"
    sed -i "s/relax_Y_minq = .*/relax_Y_minq = $result/" inlist_hb
    sed -i "s/relax_Z_minq = .*/relax_Z_minq = $result/" inlist_hb

    sed -i "s#load_model_filename = .*#load_model_filename = './mod_env/ZAHB_$1-$3-$2.mod'#" inlist_hb
    sed -i "s#save_model_filename = .*#save_model_filename = './mod_TAHB/TAHB_$1-$3-$4.mod'#" inlist_hb

    sed -i "s/initial_mass = .*/initial_mass = $1/" inlist_hb
    sed -i "s/Zbase = .*/Zbase = $2/" inlist_hb
    sed -i "s/initial_z = .*/initial_z = $2/" inlist_hb
    sed -i "s/extra_mass_retained_by_remove_H_env = .*/extra_mass_retained_by_remove_H_env = $3/" inlist_hb
    sed -i "s/new_Y = .*/new_Y = $4/" inlist_hb

    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
    do_one inlist_hb_header
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

    cp LOGS_HB/history.data data_hb/$1-$3-$4.data
    rm -f ./LOGS_HB/*
}

function modify_inlist_wd {
    echo -e "\nZbase=$new_Zbase,M=$new_mass,env_mass=$new_env_mass,Y_surf=$Y_surf"
    

    sed -i "s#load_model_filename = .*#load_model_filename = './mod_TAHB/TAHB_$1-$3-$4.mod'#" inlist_wd

    sed -i "s/initial_mass = .*/initial_mass = $1/" inlist_wd
    sed -i "s/Zbase = .*/Zbase = $2/" inlist_wd
    sed -i "s/initial_z = .*/initial_z = $2/" inlist_wd

    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
    do_one inlist_wd_header
    date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

    cp LOGS_WD/history.data data_wd/$1-$3-$4.data
    rm -f ./LOGS_WD/*
}

source "${MESA_DIR}/star/test_suite/test_suite_helpers"

new_Zbase="0.02"




for new_mass in $(seq 4.8 0.1 6.0)
do
    # modify_inlist_to_TAMS $new_mass $new_Zbase
    # modify_inlist_to_ZACHeB $new_mass $new_Zbase
    for new_env_mass in 0.0001 0.0005 0.0010 0.0020 0.0050 0.0070 0.0100 0.0200
    do
        # modify_inlist_remove_env $new_mass $new_Zbase $new_env_mass
        for Y_surf in $(seq 0.10 0.10 0.90)
        do
            # modify_inlist_hb $new_mass $new_Zbase $new_env_mass $Y_surf
            modify_inlist_wd $new_mass $new_Zbase $new_env_mass $Y_surf
        done
        # rm -f ./LOGS_env/*
    done
done
#$(seq 2.4 0.1 6.0)
# export env_mass
# env_mass2=$(echo "$new_env_mass + $new_env_mass" | bc)
# echo $env_mass2
# 0.0001 0.0005 0.0010 0.0020 0.0050 0.0070 0.0100 0.0200 0.0500 0.0700 0.1000






