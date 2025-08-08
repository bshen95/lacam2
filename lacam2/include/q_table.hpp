#pragma once
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <limits>
#include <utility>


using State  = std::size_t;
using Action = std::size_t;

/// Hash for (state, action) pairs when the table is in sparse mode
struct QPairHash {
    std::size_t operator()(const std::pair<State,Action>& p) const noexcept
    {
        return std::hash<State>{}(p.first) ^ (std::hash<Action>{}(p.second) << 1);
    }
};

class QTable
{
public:
    /// Construct a dense table.  Pass 0 for either dimension to switch to sparse mode.
    QTable(std::size_t nStates    = 0,
            std::size_t nActions   = 0,
            double      initValue  = 0.0)
    : dense_{nStates && nActions}
    , nStates_{nStates}
    , nActions_{nActions}
    {
        if (dense_)
            data_.assign(nStates_ * nActions_, initValue);
        else
            sparseInitValue_ = initValue;
    }

    /// Read Q(s,a)
    double operator()(State s, Action a) const
    {
        if (dense_) {
            return data_[idx(s,a)];
        }
        const auto it = sparse_.find({s,a});
        return it == sparse_.end() ? sparseInitValue_ : it->second;
    }

    /// Set Q(s,a)
    void set(State s, Action a, double value)
    {
        if (dense_) {
            data_[idx(s,a)] = value;
        } else {
            sparse_[{s,a}] = value;
        }
    }

    /// Incremental update: Q(s,a) ← Q(s,a) + α·δ
    //Q(s,a)←Q(s,a)+α⋅(Gt − Q(s,a))
    void update(State s, Action a, double alpha, double delta)
    {
        if (dense_) {
            data_[idx(s,a)] += alpha * delta;
        } else {
            sparse_[{s,a}] = (*this)(s,a) + alpha * delta;
        }
    }

    void dump_q_table_csv(const std::string& filename) const
    {
        std::ofstream fout(filename);
        fout << "state,action,q_value\n";
        if (dense_) {
            for (State s = 0; s < nStates_; ++s) {
                for (Action a = 0; a < nActions_; ++a) {
                    fout << s << "," << a << "," << data_[idx(s,a)] << "\n";
                }
            }
        } else {
            // For sparse, output all possible state-action pairs (including default values)
            for (State s = 0; s < nStates_; ++s) {
                for (Action a = 0; a < nActions_; ++a) {
                    fout << s << "," << a << "," << (*this)(s,a) << "\n";
                }
            }
        }
        fout.close();
    }
    
    std::vector<std::tuple<State, int, double>> get_q_value_records(State s)
    {
        std::vector<std::tuple<State, int, double>> records;
        if (dense_) {
            for (Action a = 0; a < nActions_; ++a) {
                records.emplace_back(s, a, data_[idx(s,a)]);
            }
        } else {
            for (Action a = 0; a < nActions_; ++a) {
                records.emplace_back(s, a, (*this)(s,a));
            }
        }
        return records;
    }

    void print_Q_value_at_state(State s)
    {
        if (dense_) {
            for (Action a = 0; a < nActions_; ++a) {
                std::cout << "Q(" << s << "," << a << ") = " << data_[idx(s,a)] << "\n";
            }
        } else {
            for (Action a = 0; a < nActions_; ++a) {
                std::cout << "Q(" << s << "," << a << ") = " << (*this)(s,a) << "\n";
            }
        }
    }

    /// Arg‑max action for a state (dense mode) or using `actionCountHint` in sparse mode.
    Action bestAction(State s, std::size_t actionCountHint = 0) const
    {
        if (dense_) {
            double bestVal = -std::numeric_limits<double>::infinity();
            Action bestAct = 0;
            for (Action a = 0; a < nActions_; ++a)
            {
                double q = data_[idx(s,a)];
                if (q > bestVal) { bestVal = q; bestAct = a; }
            }
            return bestAct;
        }
        double bestVal = -std::numeric_limits<double>::infinity();
        Action bestAct = 0;
        for (Action a = 0; a < actionCountHint; ++a)
        {
            double q = (*this)(s,a);
            if (q > bestVal) { bestVal = q; bestAct = a; }
        }
        return bestAct;
    }

    // Convenience accessors
    std::size_t numStates()  const { return nStates_;  }
    std::size_t numActions() const { return nActions_; }
    bool        isDense()    const { return dense_;    }

private:
    bool dense_;
    std::size_t nStates_{0}, nActions_{0};

    // Dense representation
    std::vector<double> data_;

    // Sparse representation
    std::unordered_map<std::pair<State,Action>, double, QPairHash> sparse_;
    double sparseInitValue_{0.0};

    // Helpers
    std::size_t idx(State s, Action a) const noexcept { return s * nActions_ + a; }
};
