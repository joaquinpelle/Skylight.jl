"""
    tmap(f, itr; tasks_per_thread::Int = 2, kwargs...)

Apply a given function `f` to partitions of the iterable `itr` in parallel.

# Arguments
- `f`: A function to be applied to each partition of the iterable.
- `itr`: The iterable to be partitioned and processed.
- `tasks_per_thread::Int`: An optional parameter to control the number of tasks per thread. Default value is 2.
- `kwargs...`: Additional keyword arguments that will be passed to the function `f`.

# Description
This function divides the iterable `itr` into chunks, with each chunk size determined by the given `tasks_per_thread` and the number of available threads. The specified function `f` is then applied to each chunk in parallel using the `@spawn` macro. The `@sync` macro ensures that all tasks are completed before moving on.

# Example
```julia
function print_elements(chunk)
    for i in chunk
        println(i)
    end
end
tmap(print_elements, 1:10; tasks_per_thread = 2)
"""
function tmap(f, itr, args...; tasks_per_thread::Int = 2, kwargs...)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = max(1, length(itr) ÷ (tasks_per_thread * nthreads()))
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(Iterators.partition(itr, chunk_size)) do chunk
        @spawn f(chunk, args...; kwargs...)
    end
end
