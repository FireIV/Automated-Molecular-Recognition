from rcsbapi.search import AttributeQuery


q = AttributeQuery(
    attribute="rcsb_id",
    operator="exact_match",
    value="B40",
    service="text_chem"
)
a = list(q())
print(a)
