const int max_threads = 20;

template <typename OutputType, typename InputType, typename ExecuteT> 
std::vector<OutputType> uni_run(const std::vector<InputType>& inputs, const ExecuteT& executor) {
    // assume default-constructible for simplicity
    std::vector<OutputType> outputs(inputs.size());
    
    for (size_t i = 0; i < inputs.size(); ++i) {
        outputs[i] = executor(inputs[i]);
    }
    return outputs;
}

template <typename OutputType, typename InputType, typename ExecuteT> 
std::vector<OutputType> pa_run(const std::vector<InputType>& inputs, const ExecuteT& executor) {
    // assume default-constructible for simplicity
    std::vector<OutputType> outputs(inputs.size());

    std::vector<int> tasks(inputs.size());
    std::iota(ALL(tasks), 0);
    std::reverse(ALL(tasks));

    std::mutex get_task_mutex;

    auto get_task = [&]() -> optional<int> {
        std::lock_guard<std::mutex> guard(get_task_mutex);
        if (tasks.empty()) {
            return {};
        }
        int result = tasks.back();
        tasks.pop_back();
        return {result};
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < max_threads; ++i) {
        threads.push_back(std::thread([&]() {
            while (true) {
                auto task = get_task();
                if (!task.has_value()) {
                    return;
                }

                outputs[task.value()] = executor(inputs[task.value()]);
                std::cerr << ".";
            }
        }));
    }

    for (int i = 0; i < max_threads; ++i) {
        threads[i].join();
    }

    return outputs;
}

struct Output {
    int result{};

    void print() const {
        std::cout << result << "\n";
    }
};


struct Input {
    static Input read() {
        Input result;
        return result;
    }
};

// usage: pa_run(inputs, f), where f: const Input& -> Output.
