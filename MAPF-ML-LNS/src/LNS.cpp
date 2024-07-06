#include "LNS.h"
#include <queue>

LNS::LNS(const Instance& instance, double time_limit, string init_algo_name, string replan_algo_name, string destory_name,
         int neighbor_size, int num_of_iterations, int screen, PIBTPPS_option pipp_option) :
         instance(instance), time_limit(time_limit), init_algo_name(std::move(init_algo_name)),
         replan_algo_name(replan_algo_name), neighbor_size(neighbor_size), num_of_iterations(num_of_iterations),
         screen(screen), path_table(instance.map_size),pipp_option(pipp_option), replan_time_limit(time_limit / 100)
{
    start_time = Time::now();
    //replan_time_limit=2;
    if (destory_name == "Adaptive")
    {
        ALNS = true;
        destroy_weights.assign(DESTORY_COUNT * num_neighbor_sizes, 1); 
        destroy_weights[DESTORY_COUNT * num_neighbor_sizes-1]=0;
    }
    else if (destory_name == "RandomWalk")
        destroy_strategy = RANDOMWALK;
    else if (destory_name == "Intersection")
        destroy_strategy = INTERSECTION;
    else
    {
        cerr << "Destroy heuristic " << destory_name << " does not exists. " << endl;
        exit(-1);
    }

    int N = instance.getDefaultNumberOfAgents();
    agents.reserve(N);
    for (int i = 0; i < N; i++)
        agents.emplace_back(instance, i);
    preprocessing_time = ((fsec)(Time::now() - start_time)).count();
    if (screen >= 2)
        cout << "Pre-processing time = " << preprocessing_time << " seconds." << endl;
}

void LNS::printAgentPath(string file_name)const
{
    freopen(file_name.c_str(),"a",stdout);
    int N = instance.getDefaultNumberOfAgents();
    printf("%d\n",N);
    vector<int> p;
    for(int i=0;i<N;i++)
    {
        p.clear();
        p.resize(instance.tot_cell,0);
        //cerr<<instance.tot_cell<<endl;
        //cerr<<agents[i].path.size()<<endl;
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            int row=loc/instance.num_of_rows;
            int col=loc%instance.num_of_rows;
            //cerr<<loc<<" "<<instance.cell_index[row][col]<<endl;

            p[instance.cell_index[row][col]]=1;
        }
        printf("[%d",p[0]);
        for(int j=1;j<instance.tot_cell;j++)printf(",%d",p[j]);
        printf("]\n");
    }
    fclose(stdout);
}

vector<double> LNS::getAgentFeature(int x,int y,int k,int t)const
{
    vector<double> feat;
    feat.clear();
    if(k<0)feat.resize(6,0);
    else
    {
        int s_loc=agents[k].path_planner.start_location;
        int g_loc=agents[k].path_planner.goal_location;
        double sx=s_loc/instance.num_of_rows,sy=s_loc%instance.num_of_rows;
        double gx=g_loc/instance.num_of_rows,gy=g_loc%instance.num_of_rows;
        feat.push_back((x-sx)/instance.num_of_rows);
        feat.push_back((y-sy)/instance.num_of_cols);
        feat.push_back((gx-x)/instance.num_of_rows);
        feat.push_back((gy-y)/instance.num_of_cols);
        int n_loc=t<agents[k].path.size()-1?agents[k].path[t+1].location:-1;
        if(n_loc<0)feat.push_back(-1),feat.push_back(-1);
        else
        {
            int nx=n_loc/instance.num_of_rows,ny=n_loc%instance.num_of_rows;
            feat.push_back(nx-x),feat.push_back(ny-y);
        }
        
    }
    return feat;
}

void printlist(vector<double> &a)
{
	printf("[ %.4lf", a[0]);
	for(int i=1;i<int(a.size());i++)printf(",%.4lf",a[i]);
	printf("]\n");
}

void printlist(vector<int> &a)
{
	printf("[ %d", a[0]);
	for(int i=1;i<int(a.size());i++)printf(",%d",a[i]);
	printf("]\n");
}

void LNS::printAgentFeatures(string file_name)const
{
    //cerr<<"here"<<endl;
    int N = instance.getDefaultNumberOfAgents();
    vector<int> p;
    p.clear();
    p.resize((instance.num_of_rows+1)*(instance.num_of_cols+1),0);
    for(int i=0;i<N;i++)
    {    
        //cerr<<instance.tot_cell<<endl;
        //cerr<<agents[i].path.size()<<endl;
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            //int row=loc/instance.num_of_rows;
            //int col=loc%instance.num_of_rows;
            //cerr<<loc<<" "<<instance.cell_index[row][col]<<endl;
            p[loc]++;
        }
    }
    vector<vector<double> > features;
    features.clear();
    for (int i = 0; i < agents.size(); i++)
    {
        vector<double> agentFeatures;
        agentFeatures.clear();
        int delays = agents[i].getNumOfDelays();
        agentFeatures.push_back(delays);
        agentFeatures.push_back(agents[i].getIShortestPath());
        agentFeatures.push_back(agents[i].getCost());
        agentFeatures.push_back(delays*1./agents[i].getIShortestPath());
        int S=0,Max=0,Min=1e9;
        vector<int> degList;
        degList.resize(6,0);
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            S+=p[loc];
            Max=max(Max,p[loc]);
            Min=min(Min,p[loc]);
            degList[instance.getDegree(loc)]++;
        }
        agentFeatures.push_back(Min);
        agentFeatures.push_back(Max);
        agentFeatures.push_back(S);
        agentFeatures.push_back(S/double(agents[i].path.size()));
        agentFeatures.push_back(degList[1]);
        agentFeatures.push_back(degList[2]);
        agentFeatures.push_back(degList[3]);
        agentFeatures.push_back(degList[4]);
        agentFeatures.push_back(agents[i].getStartRow());
        agentFeatures.push_back(agents[i].getStartCol());
        agentFeatures.push_back(agents[i].getGoalRow());
        agentFeatures.push_back(agents[i].getGoalCol());
        agentFeatures.push_back(instance.getDegree(agents[i].getGoal()));
        features.push_back(agentFeatures);
    }
    vector<double> norm;norm.clear();
    for(int i=0;i<int(features[0].size());i++)
    {
        double x=0;
        for (auto f:features)
            x=max(x,abs(f[i]));
        norm.push_back(x);
    }
    freopen(file_name.c_str(),"a",stdout);
    printf("%d\n",agents.size());
    for(auto f:features)
    {
                for(int i=0;i<int(f.size());i++)
            if(norm[i]>0)f[i]/=norm[i];
        printlist(f);
    }
    fclose(stdout);
    //cerr<<"here"<<endl;
}

void LNS::printNodeFeature(string file_name)const
{
    freopen(file_name.c_str(),"a",stdout);
    int tp=-1;
    printf("%d\n",instance.tot_cell);
    for(int i=0;i<instance.num_of_rows;i++)
        for(int j=0;j<instance.num_of_cols;j++)
        {
            tp++;
            vector<double> feat;feat.clear();
            if(instance.cell_index[i][j]<0)continue;
            //cerr<<path_table.makespan<<endl;
            for(int k=0;k<int(path_table.table[tp].size());k++)
            {
                //cerr<<i<<" "<<j<<" "<<k<<" "<<path_table.table[tp][k]<<endl;
                vector<double> agent_feature=getAgentFeature(i,j,path_table.table[tp][k],k);
                for(auto x:agent_feature)feat.push_back(x);
            }
            if(feat.size()==0)
            {
                printf("[]\n");
                continue;
            }
            printf("[%.5lf",feat[0]);
            for(int k=1;k<int(feat.size());k++)printf(",%.5lf",feat[k]);
            printf("]\n");
        }
    
    fclose(stdout);
}

void LNS::addEdgetoDG(string file_name)
{
    //cerr<<"here2"<<endl;

    int N = instance.getDefaultNumberOfAgents();
    freopen((file_name+"GraphDG").c_str(),"r",stdin);
    timestamp++;
    int _,__;
    vector<int> source,target,weight;
    source.clear(),target.clear(),weight.clear();
    int new_edge=0;
    scanf("%d",&_);for(;_;_--){scanf("%d",&__);source.push_back(__);}
    scanf("%d",&_);for(;_;_--){scanf("%d",&__);target.push_back(__);}
    scanf("%d",&_);for(;_;_--){scanf("%d",&__);weight.push_back(__);}
    fclose(stdin);
    for(int i=0;i<(int)source.size();i++)AG[source[i]][target[i]]=timestamp;
    //cerr<<"here2"<<endl;
    
    for (int a=0;a<N;a++)
    {
        set<int> neighbors_set;
        neighbors_set.clear();
        neighbors_set.insert(a);
        for(_=0;_<6;_++)randomWalk(a, agents[a].path[0].location, 0, neighbors_set, N*N, (int) agents[a].path.size() - 1);
        set<int>::iterator it;
        for (it=neighbors_set.begin(); it!=neighbors_set.end(); ++it)
        {
            if(AG[a][*it]==timestamp||a==*it)continue;
            AG[a][*it]=timestamp;
            source.push_back(a);
            target.push_back(*it);
            weight.push_back(0);
            new_edge++;
        }
    }
    //cerr<<"here2"<<endl;

    freopen((file_name+"Graph").c_str(),"w",stdout);
    printlist(source);
    printlist(target);
    printlist(weight);
    fclose(stdout);
    //cerr<<new_edge<<" "<<source.size()<<" ";
}

void LNS::predictNeighbor()
{
    cerr<<"ML"<< " ";
    if (num_of_model>1)
    {
        system(("python3 learning3/predict.py --model "+model+std::to_string(selected_model)+" --feature "+feature_file).c_str());
    }
    else 
        system(("python3 learning3/predict.py --model "+model+" --feature "+feature_file).c_str());
    //cerr<<"here"<<endl;

    freopen((feature_file+"_predict").c_str(),"r",stdin);
    int len_neigh;
    scanf("%d",&len_neigh);
    //cerr<<initial_sum_of_costs<<endl;
    neighbor.agents.clear();
    for(int i=0;i<len_neigh;i++)
    {
        int x;scanf("%d",&x);
        //cerr<<x<<" ";
        neighbor.agents.push_back(x);
    }
    fclose(stdin);
}


void LNS::printAgentFeatures2(string feature_file,int id)const
{
    //cerr<<"here"<<endl;
    int N = instance.getDefaultNumberOfAgents();
    vector<int> p;
    p.clear();
    p.resize((instance.num_of_rows+1)*(instance.num_of_cols+1),0);
    for(int i=0;i<N;i++)
    {    
        //cerr<<instance.tot_cell<<endl;
        //cerr<<agents[i].path.size()<<endl;
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            //int row=loc/instance.num_of_rows;
            //int col=loc%instance.num_of_rows;
            //cerr<<loc<<" "<<instance.cell_index[row][col]<<endl;
            p[loc]++;
        }
    }
    vector<int> flag;flag.clear();
    flag.resize(N,0);
    for(auto x: neighbor.agents)flag[x]=1;
    vector<vector<double> > features;
    features.clear();
    for (int i = 0; i < agents.size(); i++)
    {
        vector<double> agentFeatures;
        agentFeatures.clear();
        int delays = agents[i].getNumOfDelays();
        agentFeatures.push_back(flag[i]);
        agentFeatures.push_back(delays);
        agentFeatures.push_back(agents[i].getIShortestPath());
        agentFeatures.push_back(agents[i].getCost());
        agentFeatures.push_back(delays*1./agents[i].getIShortestPath());
        int S=0,Max=0,Min=1e9;
        vector<int> degList;
        degList.resize(6,0);
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            S+=p[loc];
            Max=max(Max,p[loc]);
            Min=min(Min,p[loc]);
            degList[instance.getDegree(loc)]++;
        }
        agentFeatures.push_back(Min);
        agentFeatures.push_back(Max);
        agentFeatures.push_back(S);
        agentFeatures.push_back(S/double(agents[i].path.size()));
        agentFeatures.push_back(degList[1]);
        agentFeatures.push_back(degList[2]);
        agentFeatures.push_back(degList[3]);
        agentFeatures.push_back(degList[4]);
        agentFeatures.push_back(agents[i].getStartRow());
        agentFeatures.push_back(agents[i].getStartCol());
        agentFeatures.push_back(agents[i].getGoalRow());
        agentFeatures.push_back(agents[i].getGoalCol());
        agentFeatures.push_back(instance.getDegree(agents[i].getGoal()));
        features.push_back(agentFeatures);
    }
    vector<double> norm;norm.clear();
    for(int i=0;i<int(features[0].size());i++)
    {
        double x=0;
        for (auto f:features)
            x=max(x,abs(f[i]));
        norm.push_back(x);
    }
    string file_name=feature_file;
    if(id>9)file_name+=char('0'+id/10);
    file_name+=char('0'+id%10);
    freopen(file_name.c_str(),"w",stdout);
    printf("%d\n",agents.size());
    for(auto f:features)
    {
        for(int i=0;i<int(f.size());i++)
            if(norm[i]>0)f[i]/=norm[i];
        printlist(f);
    }
    //printf("[%d,%d,%d]\n",neighbor.old_sum_of_costs-neighbor.sum_of_costs,neighbor.old_sum_of_costs,neighbor.sum_of_costs);
    fclose(stdout);
    //cerr<<"here"<<endl;
}

void LNS::getAgentFeatures2()
{
    //cerr<<"here"<<endl;
    int N = instance.getDefaultNumberOfAgents();
    vector<int> p;
    p.clear();
    p.resize((instance.num_of_rows+1)*(instance.num_of_cols+1),0);
    for(int i=0;i<N;i++)
    {    
        //cerr<<instance.tot_cell<<endl;
        //cerr<<agents[i].path.size()<<endl;
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            //int row=loc/instance.num_of_rows;
            //int col=loc%instance.num_of_rows;
            //cerr<<loc<<" "<<instance.cell_index[row][col]<<endl;
            p[loc]++;
        }
    }
    vector<int> flag;flag.clear();
    flag.resize(N,0);
    for(auto x: neighbor.agents)flag[x]=1;
    //vector<vector<double> > features;
    Features.clear();
    for (int i = 0; i < agents.size(); i++)
    {
        vector<double> agentFeatures;
        agentFeatures.clear();
        int delays = agents[i].getNumOfDelays();
        agentFeatures.push_back(flag[i]);
        agentFeatures.push_back(delays);
        agentFeatures.push_back(agents[i].getIShortestPath());
        agentFeatures.push_back(agents[i].getCost());
        agentFeatures.push_back(delays*1./agents[i].getIShortestPath());
        int S=0,Max=0,Min=1e9;
        vector<int> degList;
        degList.resize(6,0);
        for(int j=0;j<int(agents[i].path.size());j++)
        {
            int loc=agents[i].path[j].location;
            S+=p[loc];
            Max=max(Max,p[loc]);
            Min=min(Min,p[loc]);
            degList[instance.getDegree(loc)]++;
        }
        agentFeatures.push_back(Min);
        agentFeatures.push_back(Max);
        agentFeatures.push_back(S);
        agentFeatures.push_back(S/double(agents[i].path.size()));
        agentFeatures.push_back(degList[1]);
        agentFeatures.push_back(degList[2]);
        agentFeatures.push_back(degList[3]);
        agentFeatures.push_back(degList[4]);
        agentFeatures.push_back(agents[i].getStartRow());
        agentFeatures.push_back(agents[i].getStartCol());
        agentFeatures.push_back(agents[i].getGoalRow());
        agentFeatures.push_back(agents[i].getGoalCol());
        agentFeatures.push_back(instance.getDegree(agents[i].getGoal()));
        Features.push_back(agentFeatures);
        if (i==0)
        {
            for(int j=0;j<int(agentFeatures.size());j++)
            {
                min_feat[j].clear(),max_feat[j].clear(),sum_feat[j]=agentFeatures[j];
                min_feat[j].push_back(make_pair(agentFeatures[j],i));
                max_feat[j].push_back(make_pair(-agentFeatures[j],i));
            }
        }else{
            for(int j=0;j<int(agentFeatures.size());j++)
            {
                sum_feat[j]+=agentFeatures[j];
                min_feat[j].push_back(make_pair(agentFeatures[j],i));
                max_feat[j].push_back(make_pair(-agentFeatures[j],i));
            }
        }
    }
    for (int i=0;i<int(Features[0].size());i++)
    {
        nth_element(min_feat[i].begin(),min_feat[i].begin()+ub_neighborsize+2,min_feat[i].end());
        sort(min_feat[i].begin(),min_feat[i].begin()+ub_neighborsize+2);
        nth_element(max_feat[i].begin(),max_feat[i].begin()+ub_neighborsize+2,max_feat[i].end());
        sort(max_feat[i].begin(),max_feat[i].begin()+ub_neighborsize+2);
        for(int j=0;j<ub_neighborsize+1;j++)max_feat[i][j].first*=-1;
        //cerr<<ub_neighborsize<<endl;
    }
    /*
    vector<double> norm;norm.clear();
    for(int i=0;i<int(Features[0].size());i++)
    {
        double x=0;
        for (auto f:Features)
            x=max(x,abs(f[i]));
        norm.push_back(x);
    }

    for(int j=0;j<(int)agents.size();j++)
    {
        for(int i=0;i<int(Features[j].size());i++)
        {
            if(norm[i]>0)Features[j][i]/=norm[i];
        }
        //cerr<<endl;
    }*/

}

void LNS::process_neighborFeatures(vector<double> &a)
{
    a.clear();
    /*
    vector<int> flag;flag.clear();
    flag.resize(agents.size(),0);
    for(auto x: neighbor.agents)flag[x]=1;
    for(int i=0;i<2;i++)
        for(int j=0;j<(int)Features[0].size();j++)
        {
            feat_part[i][j][0]=1e9,feat_part[i][j][1]=-1e9;
            feat_part[i][j][2]=feat_part[i][j][3]=0;
        }
    for(int i=0;i<(int)agents.size();i++)
    {
        int k=flag[i];
        for(int j=1;j<int(Features[i].size());j++)
        {
            feat_part[k][j][0]=min(feat_part[k][j][0],Features[i][j]);
            feat_part[k][j][1]=max(feat_part[k][j][1],Features[i][j]);
            feat_part[k][j][2]+=Features[i][j];
            feat_part[k][j][3]++;
        }
    }*/
    a.push_back(neighbor.agents.size());
    ++timestamp;
    for(auto x: neighbor.agents)AG[x][0]=timestamp;
    for(int i=1;i<int(Features[0].size());i++)
    {
        for(int j=0;;j++)
            if(AG[min_feat[i][j].second][0]!=timestamp) 
                { a.push_back(min_feat[i][j].first); break;}
        for(int j=0;;j++)
            if(AG[max_feat[i][j].second][0]!=timestamp) 
                { a.push_back(max_feat[i][j].first); break;}
        double ss=0;
        for(auto x:neighbor.agents)ss+=Features[x][i];
        a.push_back(sum_feat[i]-ss);
        a.push_back((sum_feat[i]-ss)/(agents.size()-neighbor.agents.size()));
    }
    for(int i=1;i<int(Features[0].size());i++)
    {
        double m1=1e9,m2=-1e9;
        double ss=0;
        for(auto x:neighbor.agents)
        {
            m1=min(m1,Features[x][i]);
            m2=max(m2,Features[x][i]);
            ss+=Features[x][i];
        }
        a.push_back(m1);
        a.push_back(m2);
        a.push_back(ss);
        a.push_back(ss/neighbor.agents.size());
    }
    //for(auto x:a)cerr<<x<<" ";
    //cerr<<endl;
/*
    //cerr<<neighbor.agents.size()<<" ";
    for(int i=0;i<2;i++)
        for(int j=1;j<int(Features[0].size());j++)
        {
            cerr<<feat_part[i][j][0]<<" "<<feat_part[i][j][1]<<" "<<feat_part[i][j][2]<<" "<<feat_part[i][j][2]/feat_part[i][j][3]<<" ";
            a.push_back(feat_part[i][j][0]);
            a.push_back(feat_part[i][j][1]);
            a.push_back(feat_part[i][j][2]);
            a.push_back(feat_part[i][j][2]/feat_part[i][j][3]);
        }
    cerr<<endl;
    exit(0);*/
    //cerr<<"procend"<<endl;

}

void LNS::normalize_neighborFeatures(int num_of_samples)
{
    //freopen(feature_file.c_str(),"w",stdout);

    vector<double> norm;norm.clear();
    for(int i=0;i<int(neighbor_feat[0].size());i++)
    {
        double x=0;
        for(int j=0;j<num_of_samples;j++)
            x=max(x,abs(neighbor_feat[j][i]));
        if (abs(x)<1e-8)x=1;
        norm.push_back(x);
        //cout<<x<<" ";
    }
    //cout<<endl;
    for(int j=0;j<num_of_samples;j++)
    {
        //cout<<"h ";
        for(int i=0;i<int(neighbor_feat[j].size());i++)
        {
            neighbor_feat[j][i]/=norm[i];
            //cout<<neighbor_feat[j][i]<<" ";
        }
        //cout<<endl;
    }
    //fclose(stdout);
}

bool LNS::run()
{
    // only for statistic analysis, and thus is not included in runtime
    srand(19961012);
    //srand(19961015);
    sum_of_distances = 0;
    for (const auto & agent : agents)
    {
        sum_of_distances += agent.path_planner.my_heuristic[agent.path_planner.start_location];
    }

    //instance.getRowCoordinate
    //cerr<<"here"<<endl;
    initial_solution_runtime = 0;
    bool succ = false;
    int count = 0;
    int previous_model=-1;
    int cnt_getFeat=0;
    if (initSol_file!="NONE")
    {
        //cerr<<"loading initial solution!"<<endl;
        initial_sum_of_costs=0;
        freopen(initSol_file.c_str(),"r",stdin);
        int _;scanf("%d",&_);
        for(int i=0;i<_;i++)
            path_table.deletePath(i, agents[i].path);
        for(int i=0;i<_;i++)
        {
            int len;scanf("%d",&len);
            agents[i].path.clear();
            initial_sum_of_costs+=len;
            for(int j=0;j<=len;j++)
            {
                int loc;scanf("%d",&loc);
                agents[i].path.push_back(PathEntry(loc));
            }
            path_table.insertPath(i,agents[i].path);
        }
        int len_neigh;scanf("%d",&len_neigh);
        sum_of_costs=initial_sum_of_costs;
        //cerr<<initial_sum_of_costs<<endl;
        neighbor.agents.clear();
        for(int i=0;i<len_neigh;i++)
        {
            int x;scanf("%d",&x);
            //cerr<<x<<" ";
            neighbor.agents.push_back(x);
        }
        //cerr<<"\n";
        fclose(stdin);
        neighbor.old_paths.resize(neighbor.agents.size());
        neighbor.old_sum_of_costs = 0;
        for (int i = 0; i < (int)neighbor.agents.size(); i++)
        {
            if (replan_algo_name == "PP")
                neighbor.old_paths[i] = agents[neighbor.agents[i]].path;
            path_table.deletePath(neighbor.agents[i], agents[neighbor.agents[i]].path);
            neighbor.old_sum_of_costs += agents[neighbor.agents[i]].path.size() - 1;
        }
        if (replan_algo_name == "EECBS")
            succ = runEECBS();
        else if (replan_algo_name == "CBS")
            succ = runCBS();
        else if (replan_algo_name == "PP")
            succ = runPP();
        else
        {
            cerr << "Wrong replanning strategy" << endl;
            exit(-1);
        }
        cerr<<"! ";
        sum_of_costs+=-neighbor.old_sum_of_costs+neighbor.sum_of_costs;
        //cerr<<sum_of_costs<<" "<<neighbor.old_sum_of_costs<<" "<<neighbor.sum_of_costs<<endl;
    }else
    {
        //cerr<<"finding initial solution!"<<endl;
        while (!succ && initial_solution_runtime < time_limit)
        {
            start_time = Time::now();
            succ = getInitialSolution();
            initial_solution_runtime += ((fsec)(Time::now() - start_time)).count();
            count++;
            //cerr<<count<<endl;
        }
        iteration_stats.emplace_back(neighbor.agents.size(),
                                    initial_sum_of_costs, initial_solution_runtime, init_algo_name,0);
        runtime = initial_solution_runtime;
        if (succ)
        {
            if (screen >= 1)
                cout << "Initial solution cost = " << initial_sum_of_costs << ", "
                    << "runtime = " << initial_solution_runtime << endl;
        }
        else
        {
            cerr << "Failed to find an initial solution in "
                << runtime << " seconds and  " << count << " iterations" << endl;
            return false; // terminate because no initial solution is found
        }
        //cerr<<"found init solution!"<<endl;
    }
    //cerr<<"init_sum_of_costs"<<" "<<initial_sum_of_costs<<endl;


    if (ALNS)
    {
        //if(model!="NONE"&&num_of_model>1)
        //if(testing==1)
          //  destroy_weights.assign(3, 1);
        //for(auto x:destroy_weights)
          //  cerr<<x<<" ";
        //cerr<<endl;
    }
    isML=(testing==1&&num_of_samples>1);
    //cerr<<"begin search"<<endl;
    stage=1;
    int cnt_below=0;
    int no_improved=0;
    //srand(time(0));
    while (runtime < time_limit && iteration_stats.size() <= num_of_iterations)
    {
        //cerr<<"aaa"<<endl;
        runtime =((fsec)(Time::now() - start_time)).count();
        /*if(screen >= 1)
            validateSolution();*/
        
        double best_obj=-1;
        Neighbor best_neighbor;
        //cerr<<num_of_samples<<endl;

        ///testing part
        //if(testing==1)
        {
            if(no_improved==0){
                high_resolution_clock::time_point getFeature_starttime=Time::now();
                getAgentFeatures2();cnt_getFeat++;
                sampling_time+=((fsec)(Time::now() - getFeature_starttime)).count();
            }
            
            /*for(auto x:Features)
            {
                for(auto y:x)cerr<<y<<" ";
                cerr<<endl;
            }*/
        }
        if(cnt_below>=100)
        {
            testing=0;
            num_of_samples=1;
        }
        ///////////
        vector<double> imp;double max_imp=0;
        vector<double> Eff;
        imp.clear();imp.resize(num_of_samples,0);
        Eff.clear();Eff.resize(num_of_samples,0);
        //
        //cerr<<"ss ";


        for(int _=0;_<num_of_samples;_++)
        {

            if (ALNS)
                chooseDestroyHeuristicbyALNS();
            previous_model=selected_model;
            //if (num_of_samples==1)cerr<<selected_model<<" ";
            
            /*if(num_of_samples==1&&model!="NONE")
            {
                if (selected_model>num_of_model)
                {
                    if (selected_neighbor==num_of_model+1)succ=generateNeighborByRandomWalk();
                    else succ = generateNeighborByIntersection();
                    if(!succ)continue;
                }else
                {
                    if(0)
                    {
                        instance.printGraphFeature(feature_file);
                        //node_feature
                        //cerr<<"here!!!"<<endl;
                        printNodeFeature(feature_file);
                        //cerr<<"here!!!"<<endl;
                        //agent_path
                        printAgentPath(feature_file);
                        freopen(feature_file.c_str(),"a",stdout);
                        printf("[]\n");
                        fclose(stdout);
                    }else
                    {
                        addEdgetoDG(feature_file);
                        //printAgentFeatures2(feature_file);
                        
                        printAgentFeatures(feature_file);
                        freopen(feature_file.c_str(),"a",stdout);
                        printf("[]\n");
                        fclose(stdout);
                        
                    }
                    predictNeighbor();
                }
            }else*/
            {
                if(testing==0)
                {
                    
                    //if(_>=int(num_of_samples/3.*2)) //use random sampling
                    {
                        //cerr<<"d ";
                        switch (destroy_strategy)
                        {
                            case RANDOMWALK:
                                while(!(succ = generateNeighborByRandomWalk()));
                                break;
                            case INTERSECTION:
                                while(!(succ = generateNeighborByIntersection()));
                                break;
                            case RANDOMSAMPLE:
                            {
                                neighbor_size=int((rand()%10000/10000.)*(ub_neighborsize+1-lb_neighborsize))+lb_neighborsize;
                                //cerr<<"h ";
                                neighbor.agents.clear();
                                vector<int> idx;idx.clear();
                                for(int i=0;i<(int)agents.size();i++)idx.push_back(i);
                                random_shuffle(idx.begin(),idx.end());
                                for(int i=0;i<neighbor_size;i++)
                                {
                                    neighbor.agents.push_back(idx[i]);
                                }
                                succ=1;
                                break;
                            }
                            default:
                                cerr << "Wrong neighbor generation strategy" << endl;
                                exit(-1);
                        }
                    }
                    if(!succ)
                    {
                        cerr << "Unable to generate neighborhood" << endl;
                        exit(-1);
                        continue;
                    }
                }else
                {
                    high_resolution_clock::time_point sampling_starttime=Time::now();
                    
                    /*For Deep ML
                    addEdgetoDG(feature_file);
                    */
                    
                    neighbor_set.clear();
                    int dictator=0;
                    for(int __=0;__<num_of_samples;__++)
                    {
                        if (ALNS)
                            chooseDestroyHeuristicbyALNS();
                        //cerr<<destroy_strategy<<" ";
                        if(!dictator)
                        {
                            switch (destroy_strategy)
                            {
                                case RANDOMWALK:
                                {
                                    int cnt_rw=0;
                                    while(!(succ = generateNeighborByRandomWalk())&&cnt_rw<15)cnt_rw++;
                                    if(!succ)dictator=1;
                                    break;
                                }
                                case INTERSECTION:
                                    while(!(succ = generateNeighborByIntersection()));
                                    break;
                                case RANDOMSAMPLE:
                                {
                                    neighbor_size=int((rand()%10000/10000.)*(ub_neighborsize+1-lb_neighborsize))+lb_neighborsize;
                                    neighbor.agents.clear();
                                    vector<int> idx;idx.clear();
                                    for(int i=0;i<(int)agents.size();i++)idx.push_back(i);
                                    random_shuffle(idx.begin(),idx.end());
                                    for(int i=0;i<neighbor_size;i++)
                                    {
                                        neighbor.agents.push_back(idx[i]);
                                    }
                                    succ=1;
                                    break;
                                }
                                default:
                                    cerr << "Wrong neighbor generation strategy" << endl;
                                    exit(-1);
                            }
                        }
                        if(dictator==1)
                        {
                            while(!(succ = generateNeighborByIntersection()));
                        }
                        if(!succ)
                        {
                            cerr << "Unable to generate neighborhood" << endl;
                            exit(-1);
                            continue;
                        }
                        //cerr<<neighbor.agents.size()<<" ";
                        //cerr<<destroy_strategy<<" ";

                        /* For Deep ML
                        string file_name=feature_file;
                        if(__>9)file_name+=char('0'+__/10);
                        file_name+=char('0'+__%10);
                        freopen(file_name.c_str(),"w",stdout);
                        printf("%d\n",agents.size());
                        vector<int> flag;flag.clear();
                        flag.resize(agents.size(),0);
                        for(auto x: neighbor.agents)flag[x]=1;
                        for(int x=0;x<int(Features.size());x++)
                        {
                            Features[x][0]=flag[x];
                            printlist(Features[x]);
                        }
                        printf("[%d,%d,%d]\n",destroy_strategy,__,0);
                        fclose(stdout);*/
                        neighbor_set.push_back(neighbor);
                        process_neighborFeatures(neighbor_feat[__]);
                    }
                    //cerr<<endl;
                    sampling_time+=((fsec)(Time::now() - sampling_starttime)).count();
                }
            }
            //predicting socket
            if(testing==1)
            {
                high_resolution_clock::time_point ml_forwardpass_starttime=Time::now();
                /*For Deep ML
                //system(("python3 learning3/predict.py --model "+model+std::to_string(stage)+" --feature "+feature_file+" --sample "+std::to_string(num_of_samples)).c_str());
                int sock = 0, valread; 
                struct sockaddr_in serv_addr; 
                //std::string hello = "dammit"; 
                //int t=clock();
                if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
                { 
                    printf("\n Socket creation error \n"); 
                    return -1; 
                } 
            
                serv_addr.sin_family = AF_INET; 
                serv_addr.sin_port = htons(PORT); 
                
                // Convert IPv4 and IPv6 addresses from text to binary form 
                if(inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr)<=0)  
                { 
                    printf("\nInvalid address/ Address not supported \n"); 
                    return -1; 
                } 
            
                if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
                { 
                    printf("\nConnection Failed \n"); 
                    return -1; 
                } 
                string msg=feature_file+" "+std::to_string(num_of_samples)+" "+std::to_string(stage);
                send(sock , msg.c_str() , strlen(msg.c_str()) , 0 ); 
                valread = read( sock , buffer, 1024); 
                printf("%s\n",buffer ); 
                freopen((feature_file+"_predict").c_str(),"r",stdin);
                forwardpass_time+=((fsec)(Time::now() - ml_forwardpass_starttime)).count();*/
                normalize_neighborFeatures(num_of_samples);
                for(int __=0;__<num_of_samples;__++)
                {
                    double score=0;
                    for(int i=0;i<(int)neighbor_feat[__].size();i++)
                    {
                        score+=model_weight[i]*neighbor_feat[__][i];
                        //cerr<<neighbor_feat[__][i]<<" ";
                    }
                    //cerr<<endl;
                    neighbor_set[__].score=score;
                }
                sort(neighbor_set.begin(),neighbor_set.end());
                forwardpass_time+=((fsec)(Time::now() - ml_forwardpass_starttime)).count();
                //cerr<<"forwarded"<<endl;
            }
            int len_neigh;
            //cerr<<"planning here!!!"<<endl;
            //cerr<<initial_sum_of_costs<<endl;
            high_resolution_clock::time_point replan_start_time=Time::now();
            int num_of_samples_to_try=num_of_samples;
            for(int x=0;x<num_of_samples_to_try;x++)
            {
                //cerr<<"~"<<x<<" ";
            // store the neighbor information
                if(testing==1)
                {
                    
                    neighbor.agents.clear();
                    for(auto y:neighbor_set[x].agents)
                        neighbor.agents.push_back(y);
                    /*For Deep ML
                    scanf("%d%d",&len_neigh,&selected_neighbor);
                    for(int i=0;i<len_neigh;i++)
                    {
                        int x;scanf("%d",&x);
                        //cerr<<x<<" ";
                        neighbor.agents.push_back(x);
                    }*/
                }
                //cerr<<"planning "<<len_neigh<<" "<<selected_neighbor<<endl;

                neighbor.old_paths.resize(neighbor.agents.size());
                neighbor.old_sum_of_costs = 0;

                //replan_start_time = Time::now();
                /*int tpsoc=0;
                for(int i=0;i<int(agents.size());i++)
                {
                    tpsoc+=agents[i].path.size() - 1;
                }
                cerr<<"tpsoc "<<tpsoc<<endl;*/

                for (int i = 0; i < (int)neighbor.agents.size(); i++)
                {
                    //cerr<<neighbor.agents[]
                    if (replan_algo_name == "PP")
                        neighbor.old_paths[i] = agents[neighbor.agents[i]].path;
                    path_table.deletePath(neighbor.agents[i], agents[neighbor.agents[i]].path);
                    neighbor.old_sum_of_costs += agents[neighbor.agents[i]].path.size() - 1;
                }

                if (replan_algo_name == "EECBS")
                    succ = runEECBS(num_of_samples>1&&testing==0);
                else if (replan_algo_name == "CBS")
                    succ = runCBS(num_of_samples>1&&testing==0);
                else if (replan_algo_name == "PP")
                    succ = runPP(num_of_samples>1&&testing==0);
                else
                {
                    cerr << "Wrong replanning strategy" << endl;
                    exit(-1); 
                }
                //cerr<<"planning "<<x<<endl;
                //cerr<<sum_of_costs<<" "<<neighbor.sum_of_costs <<" "<<neighbor.old_sum_of_costs<<endl;
                if (ALNS) // update destroy heuristics
                {
                    if (neighbor.old_sum_of_costs > neighbor.sum_of_costs )
                    {
                        destroy_weights[selected_neighbor] =
                                reaction_factor * (neighbor.old_sum_of_costs - neighbor.sum_of_costs) / neighbor.agents.size()
                                + (1 - reaction_factor) * destroy_weights[selected_neighbor];
                        previous_model=-1;
                        cnt_below=0;
                        no_improved=0;
                        //cerr<<neighbor.old_sum_of_costs - neighbor.sum_of_costs<<" ";
                        if(testing==1)break;
                    }
                    else
                    {
                        destroy_weights[selected_neighbor] =
                                (1 - decay_factor) * destroy_weights[selected_neighbor];
                        no_improved=1;
                        if(testing==1&&x<num_of_samples_to_try-1)
                        {
                            iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, ((fsec)(Time::now() - start_time)).count(), replan_algo_name,selected_model);
                            if(iteration_stats.size()>num_of_iterations)break;
                            cnt_below++;
                            if(cnt_below>=5)stage=2;
                        }
                    }
                }
                if (neighbor.old_sum_of_costs > neighbor.sum_of_costs &&testing==1)break;
                if(testing==0&&num_of_samples>1)
                {
                    printAgentFeatures2(feature_file,_);
                    imp[_]=neighbor.old_sum_of_costs-neighborhood_ave_cost/6;
                    Eff[_]=imp[_]/((fsec)(Time::now() - replan_start_time)).count();
                    max_imp=max(max_imp,imp[_]);
                    //cerr<<neighbor.old_sum_of_costs <<" "<< neighbor.sum_of_costs<<endl;
                    break;
                }
            }
            replan_time+=((fsec)(Time::now() - replan_start_time)).count();
            if(testing==1){
                fclose(stdin);
                break;
            }
            if(num_of_samples>1&&testing==0)
            {
                double replan_runtime = ((fsec)(Time::now() - replan_start_time)).count();
                
                double eff=double(neighbor.old_sum_of_costs - neighbor.sum_of_costs)/double(replan_runtime);
                //Eff[_]=eff;
                //cerr<<neighbor.old_sum_of_costs - neighbor.sum_of_costs<<" "<<replan_runtime<<" "<<eff<<"||";
                //double eff=double(neighbor.old_sum_of_costs - neighbor.sum_of_costs)/double(1);

                if(eff>best_obj && neighbor.old_sum_of_costs > neighbor.sum_of_costs)
                {
                    best_obj=eff;
                    best_neighbor=neighbor;
                }
                neighbor.sum_of_costs = neighbor.old_sum_of_costs;
                if (_%5==0)
                {
                    freopen("eff.txt","a",stdout);
                    printf("%d,",best_neighbor.old_sum_of_costs-best_neighbor.sum_of_costs);
                    fclose(stdout);
                }
            }
            
        }
        //cerr<<"end of loop"<<endl;
        if(testing==0&&num_of_samples>1)
        {
            for(int _=0;_<num_of_samples;_++)
            {
                string file_name=feature_file;
                if(_>9)file_name+=char('0'+_/10);
                file_name+=char('0'+_%10);
                freopen(file_name.c_str(),"a",stdout);
                //printf("[%.5lf,%.5lf,%d]\n",max_imp<1e-5?imp[_]:(imp[_]/max_imp>0.8?1:0),imp[_],1);
                printf("[%.5lf,%.5lf,%d]\n",Eff[_],imp[_],_);
                fclose(stdout);
            }
        }
        //cerr<<"begin output"<<endl;
        //cerr<<outputSol_file<<endl;
        if (num_of_samples>1&&testing==0)
        {
            freopen("eff.txt","a",stdout);
            printf("\n");
            fclose(stdout);
        }
        runtime = ((fsec)(Time::now() - start_time)).count();
        if(num_of_samples>1&&testing==0)
        {
            
            {
                if (0)
                {
                    instance.printGraphFeature(feature_file);
                    //node_feature
                    //cerr<<"here!!!"<<endl;
                    printNodeFeature(feature_file);
                    //cerr<<"here!!!"<<endl;
                    //agent_path
                    printAgentPath(feature_file);
                }else
                {
                    //addEdgetoDG(feature_file);
                    //printAgentFeatures(feature_file);

                }
                if(0==1)
                {
                    freopen(feature_file.c_str(),"a",stdout);
                    if (best_obj>-1)
                    {
                        vector<int> label;
                        label.clear();label.resize(agents.size(),0);
                        for(auto x:best_neighbor.agents)label[x]=1;
                        //print label
                        printf("[%d",label[0]);
                        for(int k=1;k<int(label.size());k++)printf(",%d",label[k]);
                        printf("]\n");
                    }else
                    {
                        printf("[]\n");
                    }
                    fclose(stdout);
                    cerr<<" "<<best_neighbor.agents.size();
                }
            }
            //cerr<<outputSol_file<<endl;
            freopen(outputSol_file.c_str(),"w",stdout);
            printf("%d\n",int(agents.size()));
            for(int i=0;i<int(agents.size());i++)
            {
                printf("%d",agents[i].path.size()-1);
                for(auto l:agents[i].path)printf(" %d",l.location);
                printf("\n");
            }
            fclose(stdout);
            break;
        }
        sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
        //cerr<<"sum_of_cost"<<" "<<sum_of_costs<<endl;
        if (screen >= 1)
            cout << "Iteration " << iteration_stats.size() << ", "
                 << "group size = " << neighbor.agents.size() << ", "
                 << "solution cost = " << sum_of_costs << ", "
                 << "remaining time = " << time_limit - runtime << endl;
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name,selected_model);
    }
    //cerr<<cnt_rep<<"\n";

    average_group_size = - iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);
    if(testing==0&&num_of_samples>1)return true;
    cerr << getSolverName() << ": Iterations = " << iteration_stats.size() << ", "
         << "solution cost = " << sum_of_costs << ", "
         << "initial solution cost = " << initial_sum_of_costs << ", "
         << "runtime = " << runtime << ", "
         << "group size = " << average_group_size << ", "
         << "failed iterations = " << num_of_failures << endl;
    return true;
}


bool LNS::getInitialSolution()
{
    neighbor.agents.resize(agents.size());
    for (int i = 0; i < (int)agents.size(); i++)
        neighbor.agents[i] = i;
    neighbor.old_sum_of_costs = MAX_COST;
    neighbor.sum_of_costs = 0;
    bool succ = false;
    if (init_algo_name == "EECBS")
        succ = runEECBS();
    else if (init_algo_name == "PP")
        succ = runPP();
    else if (init_algo_name == "PIBT")
        succ = runPIBT();
    else if (init_algo_name == "PPS")
        succ = runPPS();
    else if (init_algo_name == "winPIBT")
        succ = runWinPIBT();
    else if (init_algo_name == "CBS")
        succ = runCBS();
    else
    {
        cerr <<  "Initial MAPF solver " << init_algo_name << " does not exist!" << endl;
        exit(-1);
    }
    if (succ)
    {
        initial_sum_of_costs = neighbor.sum_of_costs;
        sum_of_costs = neighbor.sum_of_costs;
        return true;
    }
    else
    {
        return false;
    }

}

bool LNS::runEECBS(bool old_plan)
{
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(&agents[i].path_planner);
    }

    ECBS ecbs(search_engines, path_table, screen - 1);
    ecbs.setPrioritizeConflicts(true);
    ecbs.setDisjointSplitting(false);
    ecbs.setBypass(true);
    ecbs.setRectangleReasoning(true);
    ecbs.setCorridorReasoning(true);
    ecbs.setHeuristicType(heuristics_type::WDG, heuristics_type::GLOBAL);
    ecbs.setTargetReasoning(true);
    ecbs.setMutexReasoning(false);
    ecbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    ecbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    ecbs.setSavingStats(false);
    double w;
    if (iteration_stats.empty())
        w = 2; // initial run
    else
        w = 1.1; // replan
    ecbs.setHighLevelSolver(high_level_solver_type::EES, w);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime;
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = ecbs.solve(T, 0);
    if (succ && ecbs.solution_cost < neighbor.old_sum_of_costs && !old_plan) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *ecbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = ecbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = ecbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
        }
        if(!old_plan)neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        if (!succ)
            num_of_failures++;
    }
    return succ;
}

bool LNS::runCBS(bool old_plan)
{
    if (screen >= 2)
        cout << "old sum of costs = " << neighbor.old_sum_of_costs << endl;
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(&agents[i].path_planner);
    }

    CBS cbs(search_engines, path_table, screen - 1);
    cbs.setPrioritizeConflicts(true);
    cbs.setDisjointSplitting(false);
    cbs.setBypass(true);
    cbs.setRectangleReasoning(true);
    cbs.setCorridorReasoning(true);
    cbs.setHeuristicType(heuristics_type::WDG, heuristics_type::ZERO);
    cbs.setTargetReasoning(true);
    cbs.setMutexReasoning(false);
    cbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    cbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    cbs.setSavingStats(false);
    cbs.setHighLevelSolver(high_level_solver_type::ASTAR, 1);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime; // time limit
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = cbs.solve(T, 0);
    if (succ && cbs.solution_cost < neighbor.old_sum_of_costs&&!old_plan) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *cbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = cbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = cbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
        }
        if(!old_plan)neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        if (!succ)
            num_of_failures++;
    }
    return succ;
}

bool LNS::runPP(bool old_plan)
{
    high_resolution_clock::time_point replan_start_time=Time::now();
    cnt_rep++;
    tot_node_generated=0;
    int T_replan=1;neighborhood_ave_cost=0;
    if (old_plan)T_replan=6;
    auto shuffled_agents= neighbor.agents;
    auto p = shuffled_agents.begin();
    for(int __=1;__<=T_replan;__++)
    {
        shuffled_agents = neighbor.agents;
        std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
        if (screen >= 2) {
            for (auto id : shuffled_agents)
                cout << id << "(" << agents[id].path_planner.my_heuristic[agents[id].path_planner.start_location] <<
                    "->" << agents[id].path.size() - 1 << "), ";
            cout << endl;
        }
        int remaining_agents = (int)shuffled_agents.size();
        p = shuffled_agents.begin();
        neighbor.sum_of_costs = 0;
        runtime = ((fsec)(Time::now() - start_time)).count();
        double T = time_limit - runtime; // time limit
        if (success_replan_cnt>=30)
            replan_time_limit=success_replan_time/success_replan_cnt*2;
        if (!iteration_stats.empty()) // replan
            T = min(T, replan_time_limit);
        if (testing==0&&num_of_samples>10)replan_time_limit=1;
        auto time = Time::now();
        //cerr<<T<<" "<<old_plan<<":";
        //cerr<<":";
        while (p != shuffled_agents.end() && ((fsec)(Time::now() - time)).count() < T)
        {
            //cerr<<((fsec)(Time::now() - time)).count()<<" ";
            int id = *p;
            if (screen >= 3)
                cout << "Remaining agents = " << remaining_agents <<
                    ", remaining time = " << time_limit - runtime << " seconds. " << endl
                    << "Agent " << agents[id].id << endl;
            agents[id].path = agents[id].path_planner.findOptimalPath(path_table,T-((fsec)(Time::now() - time)).count() );
            tot_node_generated+=agents[id].path_planner.num_generated;
            if (agents[id].path.empty())
            {
                break;
            }
            neighbor.sum_of_costs += (int)agents[id].path.size() - 1;
            if (neighbor.sum_of_costs >= neighbor.old_sum_of_costs)
                break;
            path_table.insertPath(agents[id].id, agents[id].path);
            remaining_agents--;
            ++p;
        }
        //cerr<<"!"<<endl;
        //cerr<<neighbor.old_sum_of_costs<<","<<neighbor.sum_of_costs<<" \n";
        if(p == shuffled_agents.end() && neighbor.sum_of_costs < neighbor.old_sum_of_costs)
            neighborhood_ave_cost+=neighbor.sum_of_costs;
        else
            neighborhood_ave_cost+=neighbor.old_sum_of_costs;
        if (T_replan>1)
        {
            /*
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                //cerr<<neighbor.agents[]
                if (replan_algo_name == "PP")
                    neighbor.old_paths[i] = agents[neighbor.agents[i]].path;
                path_table.deletePath(neighbor.agents[i], agents[neighbor.agents[i]].path);
                neighbor.old_sum_of_costs += agents[neighbor.agents[i]].path.size() - 1;
            }*/

            auto p2 = shuffled_agents.begin();
            while (p2 != p)
            {
                int a = *p2;
                path_table.deletePath(agents[a].id, agents[a].path);
                ++p2;
            }
            if (!neighbor.old_paths.empty()&&__==T_replan)
            {
                p2 = neighbor.agents.begin();
                for (int i = 0; i < (int)neighbor.agents.size(); i++)
                {
                    int a = *p2;
                    agents[a].path = neighbor.old_paths[i];
                    path_table.insertPath(agents[a].id, agents[a].path);
                    ++p2;
                }
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
    }
    //replan_time+=((fsec)(Time::now() - replan_start_time)).count();
    //cerr<<((fsec)(Time::now() - replan_start_time)).count()<<"\n";
    //ave_cost/=T_replan;
    if(old_plan)neighbor.sum_of_costs=neighborhood_ave_cost/T_replan;
    //cerr<<neighbor.old_sum_of_costs<<" "<<neighbor.sum_of_costs<<" ";
    if (p == shuffled_agents.end() && neighbor.sum_of_costs < neighbor.old_sum_of_costs&&!old_plan) // accept new paths
    {
        success_replan_cnt+=1;
        success_replan_time+=((fsec)(Time::now() - replan_start_time)).count();
        return true;
    }
    else // stick to old paths
    {
        //cerr<<"! ";
        if (old_plan)return false;

        if (p != shuffled_agents.end())
            num_of_failures++;
        auto p2 = shuffled_agents.begin();
        while (p2 != p)
        {
            int a = *p2;
            path_table.deletePath(agents[a].id, agents[a].path);
            ++p2;
        }
        if (!neighbor.old_paths.empty())
        {
            p2 = neighbor.agents.begin();
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = *p2;
                agents[a].path = neighbor.old_paths[i];
                path_table.insertPath(agents[a].id, agents[a].path);
                ++p2;
            }
        }
        if(!old_plan)neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        return false;
    }
}

bool LNS::runPPS(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    std::mt19937* MT_S = new std::mt19937(0);
    PPS solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
//    solver.WarshallFloyd();
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}
bool LNS::runPIBT(){
    auto shuffled_agents = neighbor.agents;
     std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);

    // seed for solver
    std::mt19937* MT_S = new std::mt19937(0);
    PIBT solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}

bool LNS::runWinPIBT(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    std::mt19937* MT_S = new std::mt19937(0);
    winPIBT solver(&P,pipp_option.windowSize,pipp_option.winPIBTSoft,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}

MAPF LNS::preparePIBTProblem(vector<int> shuffled_agents){

    // seed for problem and graph
    std::mt19937* MT_PG = new std::mt19937(0);

//    Graph* G = new SimpleGrid(instance);
    Graph* G = new SimpleGrid(instance.getMapFile());

    std::vector<Task*> T;
    PIBT_Agents A;

    for (int i : shuffled_agents){
        assert(G->existNode(agents[i].path_planner.start_location));
        assert(G->existNode(agents[i].path_planner.goal_location));
        PIBT_Agent* a = new PIBT_Agent(G->getNode( agents[i].path_planner.start_location));

//        PIBT_Agent* a = new PIBT_Agent(G->getNode( agents[i].path_planner.start_location));
        A.push_back(a);
        Task* tau = new Task(G->getNode( agents[i].path_planner.goal_location));


        T.push_back(tau);
        if(screen>=5){
            cout<<"Agent "<<i<<" start: " <<a->getNode()->getPos()<<" goal: "<<tau->getG().front()->getPos()<<endl;
        }
    }

    return MAPF(G, A, T, MT_PG);

}

void LNS::updatePIBTResult(const PIBT_Agents& A,vector<int> shuffled_agents){
    int soc = 0;
    for (int i=0; i<A.size();i++){
        int a_id = shuffled_agents[i];

        agents[a_id].path.resize(A[i]->getHist().size());
        int last_goal_visit = 0;
        if(screen>=2)
            std::cout<<A[i]->logStr()<<std::endl;
        for (int n_index = 0; n_index < A[i]->getHist().size(); n_index++){
            auto n = A[i]->getHist()[n_index];
            agents[a_id].path[n_index] = PathEntry(n->v->getId());

            //record the last time agent reach the goal from a non-goal vertex.
            if(agents[a_id].path_planner.goal_location == n->v->getId() 
                && n_index - 1>=0
                && agents[a_id].path_planner.goal_location !=  agents[a_id].path[n_index - 1].location)
                last_goal_visit = n_index;

        }
        //resize to last goal visit time
        agents[a_id].path.resize(last_goal_visit + 1);
        if(screen>=2)
            std::cout<<" Length: "<< agents[a_id].path.size() <<std::endl;
        if(screen>=5){
            cout <<"Agent "<<a_id<<":";
            for (auto loc : agents[a_id].path){
                cout <<loc.location<<",";
            }
            cout<<endl;
        }
        path_table.insertPath(agents[a_id].id, agents[a_id].path);
        soc += agents[a_id].path.size()-1;
    }

    neighbor.sum_of_costs =soc;
}

void LNS::chooseDestroyHeuristicbyALNS()
{
    double sum = 0;
    for (const auto& h : destroy_weights)
        sum += h;
    if (screen >= 2)
    {
        cout << "destroy weights = ";
        for (const auto& h : destroy_weights)
            cout << h / sum << ",";
    }
    double r = (double) rand() / RAND_MAX;
    double threshold = destroy_weights[0];
    selected_neighbor = 0;
    while (threshold < r * sum)
    {
        selected_neighbor++;
        threshold += destroy_weights[selected_neighbor];
    }

    switch (selected_neighbor / num_neighbor_sizes)
    {
        case 0 : destroy_strategy = RANDOMWALK; break;
        case 1 : destroy_strategy = INTERSECTION; break;
        case 2 : destroy_strategy = RANDOMSAMPLE; break;
        default : if(model=="NONE"){ cerr << "ERROR" << endl; exit(-1);}
    }
    selected_model=selected_neighbor+1;
    // neighbor_size = (int) pow(2, selected_neighbor % num_neighbor_sizes + 1);
}

bool LNS::generateNeighborByIntersection(bool temporal)
{
    if (intersections.empty())
    {
        for (int i = 0; i < instance.map_size; i++)
        {
            if (!instance.isObstacle(i) && instance.getDegree(i) > 2)
                intersections.push_back(i);
        }
    }

    set<int> neighbors_set;
    /*int location = -1;
    auto intersection_copy = intersections;
    while (neighbors_set.size() <= 1 && !intersection_copy.empty())
    {
        neighbors_set.clear();
        auto pt = intersection_copy.begin();
        std::advance(pt, rand() % intersection_copy.size());
        location = *pt;
        if (temporal)
        {
            path_table.get_agents(neighbors_set, neighbor_size, location);
        }
        else
        {
            path_table.get_agents(neighbors_set, location);
        }
        intersection_copy.erase(pt);
    }
    if (neighbors_set.size() <= 1)
        return false;*/
    neighbor_size=int((rand()%10000/10000.)*(ub_neighborsize+1-lb_neighborsize))+lb_neighborsize;
    //cerr<<neighbor_size<<" ";
    auto pt = intersections.begin();
    std::advance(pt, rand() % intersections.size());
    int location = *pt;
    path_table.get_agents(neighbors_set, neighbor_size, location);
    if (neighbors_set.size() < neighbor_size)
    {
        set<int> closed;
        closed.insert(location);
        std::queue<int> open;
        open.push(location);
        while (!open.empty() && (int) neighbors_set.size() < neighbor_size)
        {
            int curr = open.front();
            open.pop();
            for (auto next : instance.getNeighbors(curr))
            {
                if (closed.count(next) > 0)
                    continue;
                open.push(next);
                closed.insert(next);
                if (instance.getDegree(next) >= 3)
                {
                    path_table.get_agents(neighbors_set, neighbor_size, next);
                    if ((int) neighbors_set.size() == neighbor_size)
                        break;
                }
            }
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (neighbor.agents.size() > neighbor_size)
    {
        std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
        neighbor.agents.resize(neighbor_size);
    }
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by intersection " << location << endl;
    return true;
}

bool LNS::generateNeighborByRandomWalk()
{
    neighbor_size=int((rand()%10000/10000.)*(ub_neighborsize+1-lb_neighborsize))+lb_neighborsize;
    //cerr<<neighbor_size<<" ";
    if (neighbor_size >= (int)agents.size())
    {
        neighbor.agents.resize(agents.size());
        for (int i = 0; i < (int)agents.size(); i++)
            neighbor.agents[i] = i;
        return true;
    }

    int a = findMostDelayedAgent();
    if (a < 0)
        return false;
    int aa=a;
    set<int> neighbors_set;
    neighbors_set.insert(a);
    randomWalk(a, agents[a].path[0].location, 0, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
    int count = 0;
    while (neighbors_set.size() < neighbor_size && count < 10)
    {
        int t = rand() % agents[a].path.size();
        randomWalk(a, agents[a].path[t].location, t, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
        count++;
        // select the next agent randomly
        int idx = rand() % neighbors_set.size();
        int i = 0;
        for (auto n : neighbors_set)
        {
            if (i == idx)
            {
                a = i;
                break;
            }
            i++;
        }
    }
    //cerr<<count<<" "<<neighbors_set.size()<<" ";
    //if(neighbors_set.size() < lb_neighborsize)tabu_list.insert(aa);
    //tabu_list.insert(aa);
    if (neighbors_set.size() < 2)
    {
        //tabu_list.insert(a);
        return false;
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by random walks of agent " << a
             << "(" << agents[a].path_planner.my_heuristic[agents[a].path_planner.start_location]
             << "->" << agents[a].path.size() - 1 << ")" << endl;

    return true;
}

int LNS::findMostDelayedAgent()
{
    int a = -1;
    int max_delays = -1;
    for (int i = 0; i < agents.size(); i++)
    {
        if (tabu_list.find(i) != tabu_list.end())
            continue;
        int delays = agents[i].getNumOfDelays();
        if (max_delays < delays)
        {
            a = i;
            max_delays = delays;
        }
    }
    if (max_delays == 0)
    {
        tabu_list.clear();
        return -1;
    }
    tabu_list.insert(a);
    //if(num_of_samples==1)tabu_list.insert(a);
    if (tabu_list.size() == agents.size())
        tabu_list.clear();
    return a;
}

int LNS::findRandomAgent() const
{
    int a = 0;
    int pt = rand() % (sum_of_costs - sum_of_distances) + 1;
    int sum = 0;
    for (; a < (int) agents.size(); a++)
    {
        sum += agents[a].getNumOfDelays();
        if (sum >= pt)
            break;
    }
    assert(sum >= pt);
    return a;
}

// a random walk with path that is shorter than upperbound and has conflicting with neighbor_size agents
void LNS::randomWalk(int agent_id, int start_location, int start_timestep,
                     set<int>& conflicting_agents, int neighbor_size, int upperbound)
{
    int loc = start_location;
    for (int t = start_timestep; t < upperbound; t++)
    {
        auto next_locs = instance.getNeighbors(loc);
        next_locs.push_back(loc);
        while (!next_locs.empty())
        {
            int step = rand() % next_locs.size();
            auto it = next_locs.begin();
            advance(it, step);
            int next_h_val = agents[agent_id].path_planner.my_heuristic[*it];
            if (t + 1 + next_h_val < upperbound) // move to this location
            {
                path_table.getConflictingAgents(agent_id, conflicting_agents, loc, *it, t + 1);
                loc = *it;
                break;
            }
            next_locs.erase(it);
        }
        if (next_locs.empty() || conflicting_agents.size() >= neighbor_size)
            break;
    }
}


void LNS::validateSolution() const
{
    int sum = 0;
    for (const auto& a1_ : agents)
    {
        if (a1_.path.empty())
        {
            cerr << "No solution for agent " << a1_.id << endl;
            exit(-1);
        }
        else if (a1_.path_planner.start_location != a1_.path.front().location)
        {
            cerr << "The path of agent " << a1_.id << " starts from location " << a1_.path.front().location
                << ", which is different from its start location " << a1_.path_planner.start_location << endl;
            exit(-1);
        }
        else if (a1_.path_planner.goal_location != a1_.path.back().location)
        {
            cerr << "The path of agent " << a1_.id << " ends at location " << a1_.path.back().location
                 << ", which is different from its goal location " << a1_.path_planner.goal_location << endl;
            exit(-1);
        }
        for (int t = 1; t < (int) a1_.path.size(); t++ )
        {
            if (!instance.validMove(a1_.path[t - 1].location, a1_.path[t].location))
            {
                cerr << "The path of agent " << a1_.id << " jump from "
                     << a1_.path[t - 1].location << " to " << a1_.path[t].location
                     << " between timesteps " << t - 1 << " and " << t << endl;
                exit(-1);
            }
        }
        sum += (int) a1_.path.size() - 1;
        for (const auto& a2_: agents)
        {
            if (a1_.id >= a2_.id || a1_.path.empty())
                continue;
            const auto a1 = a1_.path.size() <= a2_.path.size()? a1_ : a2_;
            const auto a2 = a1_.path.size() <= a2_.path.size()? a2_ : a1_;
            int t = 1;
            for (; t < (int) a1.path.size(); t++)
            {
                if (a1.path[t].location == a2.path[t].location) // vertex conflict
                {
                    cerr << "Find a vertex conflict between agents " << a1.id << " and " << a2.id <<
                            " at location " << a1.path[t].location << " at timestep " << t << endl;
                    exit(-1);
                }
                else if (a1.path[t].location == a2.path[t - 1].location &&
                        a1.path[t - 1].location == a2.path[t].location) // edge conflict
                {
                    cerr << "Find an edge conflict between agents " << a1.id << " and " << a2.id <<
                         " at edge (" << a1.path[t - 1].location << "," << a1.path[t].location <<
                         ") at timestep " << t << endl;
                    exit(-1);
                }
            }
            int target = a1.path.back().location;
            for (; t < (int) a2.path.size(); t++)
            {
                if (a2.path[t].location == target)  // target conflict
                {
                    cerr << "Find a target conflict where agent " << a2.id << " traverses agent " << a1.id <<
                         "'s target location " << target << " at timestep " << t << endl;
                    exit(-1);
                }
            }
        }
    }
    if (sum_of_costs != sum)
    {
        cerr << "The computed sum of costs " << sum_of_costs <<
             " is different from the sum of the paths in the solution " << sum << endl;
        exit(-1);
    }
}

void LNS::writeIterStatsToFile(string file_name) const
{
    std::ofstream output;
    output.open(file_name);
    // header
    output << "num of agents," <<
           "sum of costs," <<
           "runtime," <<
           "cost lowerbound," <<
           "MAPF algorithm," <<"Model"<< endl;

    for (const auto &data : iteration_stats)
    {
        output << data.num_of_agents << "," <<
               data.sum_of_costs << "," <<
               data.runtime << "," <<
               max(sum_of_costs_lowerbound, sum_of_distances) << "," <<
               data.algorithm << ","<<data.model<<endl;
    }
    output.close();
}

void LNS::writeResultToFile(string file_name) const
{
    std::ifstream infile(file_name);
    bool exist = infile.good();
    infile.close();
    if (!exist)
    {
        ofstream addHeads(file_name);
        addHeads << "runtime,solution cost,initial solution cost,min f value,root g value," <<
                 "iterations," <<
                 "group size," <<
                 "runtime of initial solution," <<
                 "preprocessing runtime,replan time,sampling time,forward pass time,solver name,instance name" << endl;
        addHeads.close();
    }
    ofstream stats(file_name, std::ios::app);
    stats << runtime << "," << sum_of_costs << "," << initial_sum_of_costs << "," <<
            max(sum_of_distances, sum_of_costs_lowerbound) << "," << sum_of_distances << "," <<
            iteration_stats.size() << "," << average_group_size << "," << initial_solution_runtime << "," <<
            preprocessing_time << ","<<replan_time<<","<<sampling_time<<","<<forwardpass_time<<"," << getSolverName() << "," << instance.getInstanceName() << endl;
    stats.close();
}
/*
bool LNS::generateNeighborByStart()
{
    if (start_locations.empty())
    {
        for (int i = 0; i < (int)al.agents_all.size(); i++)
        {
            auto start = ml.linearize_coordinate(al.agents_all[i].initial_location);
            start_locations[start].push_back(i);
        }
        auto it = start_locations.begin();
        while(it != start_locations.end()) // delete start locations that have only one agent
        {
            if (it->second.size() == 1)
                it = start_locations.erase(it);
            else
                ++it;
        }
    }
    if (start_locations.empty())
        return false;
    auto step = rand() % start_locations.size();
    auto it = start_locations.begin();
    advance(it, step);
    neighbors.assign(it->second.begin(), it->second.end());
    if (neighbors.size() > max_group_size ||
        (replan_strategy == 0 && neighbors.size() > group_size)) // resize the group for CBS
    {
        sortNeighborsRandomly();
        if (replan_strategy == 0)
            neighbors.resize(group_size);
        else
            neighbors.resize(max_group_size);
    }
    if (options1.debug)
        cout << "Generate " << neighbors.size() << " neighbors by start location " << it->first << endl;
    return true;
}

void LNS::sortNeighborsByRegrets()
{
    quickSort(neighbors, 0, neighbors.size() - 1, true);
    if (options1.debug) {
        for (auto agent : neighbors) {
            cout << agent << "(" << al.agents_all[agent].distance_to_goal << "->" << al.paths_all[agent].size() - 1
                 << "), ";
        }
        cout << endl;
    }
}

void LNS::sortNeighborsByStrategy()
{
    if (agent_priority_strategy == 5)
    {
        // decide the agent priority for agents at the same start location
        start_locations.clear(); // map the agents to their start locations
        for (auto i : neighbors)
            start_locations[ml.linearize_coordinate(al.agents_all[i].initial_location)].push_back(i);
        for (auto& agents : start_locations)
        {
            vector<int> agents_vec(agents.second.begin(), agents.second.end());
            quickSort(agents_vec, 0, agents_vec.size() - 1, false);
            for (int i = 0; i < (int)agents.second.size(); i++)
            {
                al.agents_all[agents_vec[i]].priority = i;
            }
        }
    }

    // sort the agents
    if (agent_priority_strategy != 0)
        quickSort(neighbors, 0, (int)neighbors.size() - 1, false);
}


void LNS::quickSort(vector<int>& agent_order, int low, int high, bool regret)
{
    if (low >= high)
        return;
    int pivot = agent_order[high];    // pivot
    int i = low;  // Index of smaller element
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or equal to pivot
        if ((regret && compareByRegrets(agent_order[j], pivot)) ||
            al.compareAgent(al.agents_all[agent_order[j]], al.agents_all[pivot], agent_priority_strategy))
        {
            std::swap(agent_order[i], agent_order[j]);
            i++;    // increment index of smaller element
        }
    }
    std::swap(agent_order[i], agent_order[high]);

    quickSort(agent_order, low, i - 1, regret);  // Before i
    quickSort(agent_order, i + 1, high, regret); // After i
}



*/