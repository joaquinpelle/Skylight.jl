export nosave, with_kw_nosave

"""
    @nosave

The `@nosave` macro is used to mark fields in a struct that should not be saved when the struct 
is serialized. It does this by wrapping the field's default value in a `NoSaveField`.

# Examples

```julia
@with_kw struct MyModel
    regular_field::Int = 1
    @nosave nosave::Int = 2
end

"""

macro nosave(ex::Expr)
    if ex.head != :(=)
        error("Invalid usage of @nosave. Expected assignment expression.")
    end

    field_name = ex.args[1] isa Symbol ? ex.args[1] : ex.args[1].args[1]
    field_type = ex.args[1] isa Expr ? esc(ex.args[1].args[2]) : Any
    init_expr = ex.args[2]

    return :( $field_name::NoSaveField{$field_type} = NoSaveField($init_expr) )
end

"""
    @with_kw_nosave

The `@with_kw_nosave` macro is used to define a new struct with keyword constructors and automatically 
override the `getproperty` method for that struct. This makes it possible to use `@nosave` fields 
transparently, as if they were regular fields. 

This macro is a combination of the `@with_kw` macro from the Parameters.jl package (which provides 
keyword constructors) and a `getproperty` method definition that unwraps `@nosave` fields when accessed.

# Examples

```julia
@with_kw_nosave struct MyModel
    regular_field::Int = 1
    @nosave nosave::Int = 2
end
```
In this example, `MyModel` is defined with a keyword constructor and a custom `getproperty` method. The
`nosave` can be accessed like a regular field, but it will not be saved to an hdf5 when a `MyModel` object is serialized.
"""

macro with_kw_nosave(ex::Expr)
    @capture(ex, struct name_{params__} <: supertype_ fields__ end)
    fields_new = prewalk(fields) do x
        if @capture(x, @nosave field_ = init_expr_)
            Expr(:(=), Expr(:(::), field.args[1], :(NoSaveField{$(field.args[2])})), 
                 :(NoSaveField($init_expr)))
        else
            x
        end
    end
    esc(quote
        @with_kw struct $(name){$(params...)} <: $(supertype) $(fields_new...) end
    end)
end