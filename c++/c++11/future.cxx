#include <iostream>
#include <future>
#include <vector>

enum class async_type_t {
    ASYNC_DEFAULT,
    ASYNC_IMM,
    ASYNC_DEFER,
    ASYNC_LAST
};

static void test_async()
{
    const std::string async_default_str = "async_default";
    const std::string async_imm_str     = "async_imm";
    const std::string async_defer_str   = "async_defer";

    auto get_str = [&async_default_str, &async_imm_str, &async_defer_str]
                   (async_type_t type) -> const std::string& {
        switch (type) {
        case async_type_t::ASYNC_DEFAULT:
            return async_default_str;
        case async_type_t::ASYNC_IMM:
            return async_imm_str;
        case async_type_t::ASYNC_DEFER:
            return async_defer_str;
        case async_type_t::ASYNC_LAST:
        default:
            throw std::range_error("invalid async type was specified");
        }
    };

    auto async_default = std::async([get_str](async_type_t type) -> std::thread::id {
        std::cout << "Hello from " << get_str(type) << " ("
                  << std::this_thread::get_id() << ")" << std::endl;
        return std::this_thread::get_id();
    }, async_type_t::ASYNC_DEFAULT);

    auto async_imm = std::async(std::launch::async,
                                [get_str](async_type_t type) -> std::thread::id {
        std::cout << "Hello from " << get_str(type) << " ("
                  << std::this_thread::get_id() << ")" << std::endl;
        return std::this_thread::get_id();
    }, async_type_t::ASYNC_IMM);

    auto async_defer = std::async(std::launch::deferred,
                                  [get_str](async_type_t type) -> std::thread::id {
        std::cout << "Hello from " << get_str(type) << " ("
                  << std::this_thread::get_id() << ")" << std::endl;
        return std::this_thread::get_id();
    }, async_type_t::ASYNC_DEFER);

    auto async_last = std::async(std::launch::deferred,
                                 [get_str](async_type_t type) {
        std::cout << "Hello from " << get_str(type) << " ("
                  << std::this_thread::get_id() << ")" << std::endl;
    }, async_type_t::ASYNC_LAST);

    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "Sleep ended (" << std::this_thread::get_id()
              << ")" << std::endl;

    std::thread::id async_default_thread_id = async_default.get();
    std::thread::id async_imm_thread_id     = async_imm.get();
    std::thread::id async_defer_thread_id   = async_defer.get();

    if (async_default_thread_id == std::this_thread::get_id()) {
        std::cout << "Default and main threads are the same" << std::endl;
    } else {
        std::cout << "Default and main threads aren't the same" << std::endl;
    }

    if (async_imm_thread_id == std::this_thread::get_id()) {
        throw std::range_error("Async and main threads have to be not the same");
    }

    if (async_defer_thread_id != std::this_thread::get_id()) {
        throw std::range_error("Deferred and main threads have to be the same");
    }

    try {
        async_last.wait();
    } catch (const std::range_error e) {
        std::cout << "Expected exception: \"" << e.what() << "\"" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Invalid exception: \"" << e.what() << "\"" << std::endl;
        std::rethrow_exception(std::make_exception_ptr(e));
    } catch (...) {
        std::cerr << "Unknow exception" << std::endl;
        std::rethrow_exception(std::current_exception());
    }
}

static std::mutex global_mutex;
static std::condition_variable global_cv;

static std::future<bool> run_task(int i, int &count) {
    auto handle = [&count](int i) -> bool {
        std::lock_guard<std::mutex> guard(global_mutex, std::adopt_lock);
        ++count;

        std::cout << "Hello from [" << i << "]: "
                  << std::this_thread::get_id() << std::endl;

        return (i % 2) == 0;
    };

    std::packaged_task<bool(int i)> task(handle);
    auto future = task.get_future();
    std::thread thread(std::move(task), i);

    thread.detach();
    return std::move(future);
}

static void test_task() {
    int num_threads = std::thread::hardware_concurrency();
    int count       = 0;

    std::vector<std::future<bool> > futures;

    for (int i = 0; i < num_threads; ++i) {
        auto future = run_task(i, count);
        futures.push_back(std::move(future));
    }

    int i = 0;
    for (auto &future : futures) {
        bool is_even = future.get();
        if (is_even != ((i % 2) == 0)) {
            throw std::range_error("Unmatched even/odd answers");
        }
        i++;
    }

    if (count != num_threads) {
        throw std::range_error("Unmatched count and num_threads");
    }
}

static std::future<int> run_promise(int find_value, int &num_completed) {
    auto my_promise         = std::make_shared<std::promise<int> >();
    std::future<int> future = my_promise->get_future();

    auto handle = [&num_completed, my_promise](int find_value) -> int {
        int i = std::thread::hardware_concurrency();
        do {
            if (i == find_value) {
                std::lock_guard<std::mutex> guard(global_mutex, std::adopt_lock);
                std::cout << "Completed promise from [" << find_value << "]: "
                          << std::this_thread::get_id() << std::endl;
                my_promise->set_value(i);
            }
        } while (i-- != 0);

        {
            std::lock_guard<std::mutex> guard(global_mutex, std::adopt_lock);
            std::cout << "Fully completed task from [" << find_value << "]: "
                      << std::this_thread::get_id() << std::endl;
            num_completed++;
        }

        global_cv.notify_one();
    };

    std::thread thread(handle, find_value);

    thread.detach();
    return std::move(future);
}

static void test_promise() {
    int num_threads   = std::thread::hardware_concurrency();
    int num_completed = 0;

    std::vector<std::future<int> > futures;

    for (int i = 0; i < num_threads; ++i) {
        auto future = run_promise(i, num_completed);
        futures.push_back(std::move(future));
    }

    int i = 0;
    for (auto &future : futures) {
        int res = future.get();
        if (res != i) {
            throw std::range_error("Unmatched answers");
        }
        i++;
    }

    std::unique_lock<std::mutex> lock(global_mutex);
    while (num_completed < num_threads) {
        if (global_cv.wait_for(lock, std::chrono::seconds(5)) ==
            std::cv_status::timeout) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
        }
    }

    if (num_completed != num_threads) {
        throw std::range_error("Unmatched count and num_threads");
    }
}

int main(int argc, char *argv[]) {
    std::cout << "Testing std::async" << std::endl;
    test_async();

    std::cout << "Testing std::packaged_task" << std::endl;
    test_task();

    std::cout << "Testing std::promise" << std::endl;
    test_promise();

    return 0;
}
