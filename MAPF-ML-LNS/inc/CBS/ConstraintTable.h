#pragma once
#include "common.h"
#include "CBSNode.h"
#include "PathTable.h"

class ConstraintTable
{
public:
	int length_min = 0;
	int length_max = MAX_TIMESTEP;
	int goal_location;
	int latest_timestep = 0; // No negative constraints after this timestep.
	size_t num_col;
	size_t map_size;

    const PathTable& path_table;

	int getHoldingTime() const; // the earliest timestep that the agent can hold its goal location

	// void clear(){ct.clear(); cat_small.clear(); cat_large.clear(); landmarks.clear(); length_min = 0, length_max = INT_MAX; latest_timestep = 0;}

	bool constrained(size_t loc, int t) const;
    bool constrained(size_t curr_loc, size_t next_loc, int next_t) const;
	int getNumOfConflictsForStep(size_t curr_id, size_t next_id, int next_timestep) const;
	// ConstraintTable() = default;
	ConstraintTable(const PathTable& path_table, size_t num_col, size_t map_size, int goal_location = -1) :
            path_table(path_table), goal_location(goal_location), num_col(num_col), map_size(map_size)
    {
        latest_timestep = max(latest_timestep, path_table.makespan);
    }
	ConstraintTable(const ConstraintTable& other) : path_table(other.path_table) {copy(other); }
    ~ConstraintTable() = default;

	void copy(const ConstraintTable& other);
	void init(const ConstraintTable& other) {copy(other); }
	void clear()
	{
	    ct.clear();
	    landmarks.clear();
	    cat.clear();
	}
	void build(const HLNode& node, int agent); // build the constraint table for the given agent at the give node 
	void buildCAT(int agent, const vector<Path*>& paths, size_t cat_size); // build the conflict avoidance table

	void insert2CT(size_t loc, int t_min, int t_max); // insert a vertex constraint to the constraint table
	void insert2CT(size_t from, size_t to, int t_min, int t_max); // insert an edge constraint to the constraint table

protected:
    // Constraint Table (CT)
	unordered_map<size_t, list<pair<int, int> > > ct; // location -> time range, or edge -> time range

	unordered_map<size_t, size_t> landmarks; // <timestep, location>: the agent must be at the given location at the given timestep

	void insertLandmark(size_t loc, int t); // insert a landmark, i.e., the agent has to be at the given location at the given timestep
	list<pair<int, int> > decodeBarrier(int B1, int B2, int t) const;

	inline size_t getEdgeIndex(size_t from, size_t to) const { return (1 + from) * map_size + to; }

private:
	vector<vector<bool> > cat; // conflict avoidance table

};

