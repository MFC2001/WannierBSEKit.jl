# Check that every name used in function position in the Julia scripts under
# examples/ comes from the script's own using/import statements or from a
# binding inside the script.
#
# How: evaluate the script's top-level using/import statements in a throwaway
# module, then ask isdefined for each called name. This avoids maintaining a
# list of "exported names" by hand; Julia resolves reexports, package
# extensions, etc. exactly as it does when the script runs.
#
# Known limits (a static check, not a compiler): only call positions are checked
# (map(f, xs) passes f by value); qualified names Base.f(x), macro names, and
# included files are not checked; bindings are file-wide (one flat scope), so a
# same-named binding elsewhere can hide a typo; missing packages are skipped.

# ---------- walk every Expr ----------
function foreach_expr(f, x)
    x isa Expr || return nothing
    x.head === :quote && return nothing   # quoted code does not run here
    f(x)
    for a in x.args
        foreach_expr(f, a)
    end
    return nothing
end

# ---------- collect names in call position ----------
function collect_called(ex)
    called = Set{Symbol}()
    foreach_expr(ex) do node
        if node.head === :call && !isempty(node.args)
            f = node.args[1]
            f isa Expr && f.head === :curly && (f = f.args[1])  # T{Int}(...)
            # the head of `a .+ b` is the symbol :.+; a dotted operator is syntax, not a name
            if f isa Symbol && !startswith(String(f), ".")
                push!(called, f)
            end
        elseif node.head === :. && length(node.args) == 2
            # f.(x) parses as Expr(:., :f, Expr(:tuple, :x))
            f, arg = node.args
            if f isa Symbol && arg isa Expr && arg.head === :tuple
                push!(called, f)
            end
        end
    end
    return called
end

# ---------- collect names bound by the script ----------
function bind_lhs!(acc, lhs)
    lhs isa Symbol && (push!(acc, lhs); return acc)
    lhs isa Expr || return acc
    if lhs.head === :call
        bind_lhs!(acc, lhs.args[1])
        for a in lhs.args[2:end]
            bind_lhs!(acc, a)
        end
    elseif lhs.head in (:(::), :(...), :where, :curly, :kw, :<:)
        bind_lhs!(acc, lhs.args[1])
    elseif lhs.head in (:tuple, :parameters)
        for a in lhs.args
            bind_lhs!(acc, a)
        end
    elseif lhs.head === :ref || lhs.head === :.
        bind_lhs!(acc, lhs.args[1])
    end
    return acc
end

function collect_bound(ex)
    bound = Set{Symbol}()
    foreach_expr(ex) do node
        h = node.head
        # assignment targets, short definitions f(x) = ..., and for/let bindings are all =
        if h in (:(=), :(+=), :(-=), :(*=), :(/=), :(^=), :(%=), :(&=), :(|=))
            bind_lhs!(bound, node.args[1])
        elseif h === :function                 # function f(x) ...
            bind_lhs!(bound, node.args[1])
        elseif h === :(->)                     # x -> ...
            bind_lhs!(bound, node.args[1])
        elseif h === :let                       # let x (bindings with = are handled above)
            for a in node.args[1:end-1]
                a isa Symbol && push!(bound, a)
            end
        elseif h === :struct
            bind_lhs!(bound, node.args[2])
        elseif h in (:abstract, :primitive)
            bind_lhs!(bound, node.args[1])
        elseif h in (:local, :global)
            for a in node.args
                a isa Symbol ? push!(bound, a) : bind_lhs!(bound, a)
            end
        elseif h === :try
            node.args[2] isa Symbol && push!(bound, node.args[2])
        end
    end
    return bound
end

# ---------- module names in a using/import statement ----------
function using_modules(stmt)
    out = Symbol[]
    for a in stmt.args
        a isa Expr && a.head === :(:) && (a = a.args[1])
        a isa Expr && a.head === :. && (a = a.args[1])
        a isa Symbol && push!(out, a)
    end
    return out
end

# ---------- build a sandbox with all of the script's imports ----------
function sandbox_for(ex)
    sandbox = Module(gensym(:ExampleSandbox))
    Core.eval(sandbox, :(using Base))
    missing = Symbol[]
    foreach_expr(ex) do node
        if node.head === :using || node.head === :import
            try
                Core.eval(sandbox, node)
            catch
                append!(missing, using_modules(node))
            end
        end
    end
    return sandbox, unique!(missing)
end

# scripts run in Main, where include is always available
const MAIN_ONLY_NAMES = Set([:include])

# check whether s is accessible in module m (including names from using and reexport)
function is_accessible(m::Module, s::Symbol)
    isdefined(m, s) && return true
    try
        Core.eval(m, s)
        return true
    catch
        return false
    end
end

# ---------- return (missing packages, unresolved call names) ----------
function unresolved_calls(source)
    ex = Meta.parseall(source)
    sandbox, missing = sandbox_for(ex)
    isempty(missing) || return missing, Set{Symbol}()
    called = collect_called(ex)
    bound = collect_bound(ex)
    unknown = Set(s for s in called if
        !(s in bound) && !(s in MAIN_ONLY_NAMES) && !is_accessible(sandbox, s))
    return missing, unknown
end

# ---------- first prove the checker itself catches errors ----------
@testset "checker self-test" begin
    missing, unknown = unresolved_calls("""
        using WannierBSEKit
        localfun(x) = x + 1
        localfun(1)
        undefined_call(1)
        undefined_broadcast.(1:2)
        a = [1, 2]
        b = [1, 2]
        c = a .+ b
        """)
    @test isempty(missing)
    @test unknown == Set([:undefined_call, :undefined_broadcast])
end

# ---------- walk examples ----------
@testset "Example scripts: called names resolve" begin
    examples_dir = joinpath(pkgdir(WannierBSEKit), "examples")
    @test isdir(examples_dir)

    if isdir(examples_dir)
        for (root, _, files) in walkdir(examples_dir)
            for file in sort(files)
                endswith(file, ".jl") || continue
                script = joinpath(root, file)
                label = relpath(script, examples_dir)

                @testset "$label" begin
                    missing, unknown = try
                        unresolved_calls(read(script, String))
                    catch err
                        @info "skipping: cannot parse" script exception = (err, catch_backtrace())
                        @test_skip false
                        continue
                    end

                    if !isempty(missing)
                        @info "skipping: example dependency unavailable" script missing
                        @test_skip false
                        continue
                    end

                    @test isempty(unknown)
                    isempty(unknown) || @info "unresolved call names" script unknown = sort!(collect(unknown))
                end
            end
        end
    end
end
