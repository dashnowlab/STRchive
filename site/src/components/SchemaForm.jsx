import { cloneElement, Fragment, useMemo, useState } from "react";
import { FaArrowDown, FaArrowUp, FaPlus, FaTrash } from "react-icons/fa6";
import { LuSend } from "react-icons/lu";
import clsx from "clsx";
import { compileSchema, draft2020, extendDraft } from "json-schema-library";
import { cloneDeep, groupBy, uniq, upperFirst } from "lodash-es";
import { get, set, split } from "@sagold/json-pointer";
import Button from "@/components/Button";
import Heading from "@/components/Heading";
import Help from "@/components/Help";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TextBox from "@/components/TextBox";
import { makeList } from "@/util/format";
import Form from "./Form";
import classes from "./SchemaForm.module.css";

/** form automatically generated from json schema */
const SchemaForm = ({ schema, data: initialData }) => {
  /** compile schema */
  const compiled = useMemo(
    () => compileSchema(schema, { drafts: [draft] }),
    [schema],
  );

  /** if no initial data, generate blank values */
  initialData ??= compiled.getData(null, { extendDefaults: false });

  /** current form data state */
  const [data, setData] = useState(initialData);

  /** revalidate when data changes */
  const { errors } = useMemo(() => {
    compiled.context.data = data;
    return compiled.validate(data);
  }, [compiled, data]);

  /** all nested fields */
  const fields = useMemo(() => {
    /** compile schema into "nodes" */
    const nodes = compiled.toSchemaNodes();

    /** return path up schema tree to root */
    const toRoot = (node) => {
      const path = [];
      while (node.parent) {
        path.push(node.parent);
        node = node.parent;
      }
      return path;
    };

    /** if section undefined, get from closest parent */
    for (const node of nodes)
      node.schema.section ??= toRoot(node).map(
        (node) => node.schema.section,
      )[0];

    /** group by section */
    return Object.entries(groupBy(nodes, "schema.section"))
      .map(([section, fields]) => ({ section, fields }))
      .filter(({ section }) => section !== "undefined");
  }, [compiled]);

  return (
    <Form className={classes.form} onSubmit={console.info}>
      {fields.map(({ section, fields }) => (
        <section key={section}>
          <Heading level={2}>{section}</Heading>

          {fields.map((field, index) => {
            /** schema props */
            const {
              type,
              title,
              description,
              examples,
              pattern,
              enum: _enum,
              minimum,
              maximum,
              less_than,
              greater_than,
            } = field.schema;

            /** path to field, without root part */
            const path = field.evaluationPath.replace(/^#\/properties/, "");

            /** normalize type to array */
            const types = [type].flat();

            /** form data name */
            const name = split(field.evaluationPath).pop();

            /** help tooltip to show */
            const tooltip = [
              description,
              ...(examples?.length ? [" "] : []),
              examples?.length ? `Examples: ${examples.join(", ")}` : "",
            ]
              .filter(Boolean)
              .join("<br/>");

            /** is the field required to be set */
            const required =
              schema.required.includes(name) && types.includes("null");

            /** full label elements to show */
            const label = (
              <>
                {title ?? name ?? ""}
                {tooltip && <Help>{tooltip}</Help>}
              </>
            );

            /** get nested data value from path */
            const value = get(data, path);

            /** set nested data value from path */
            const onChange = (value) =>
              setData(set(cloneDeep(data), path, value));

            /** get validation errors associated with this field */
            const fieldErrors = errors.filter(({ data }) =>
              data.pointer.includes(path),
            );
            /** is there an error on this field */
            const isError = !!fieldErrors.length;

            /** error component to show */
            const error = isError ? (
              <span className={classes.error}>
                {fieldErrors
                  .map(({ code, data, message }) => {
                    /** prettify error message */
                    if (code === "pattern-error")
                      return (
                        <>
                          Must match pattern <code>{data.pattern}</code>
                        </>
                      );
                    if (code === "enum-error")
                      return <>Must be {makeList(data.schema.enum, "code")}</>;
                    if (code === "compare-value-error")
                      return <>{data.schema.message}</>;
                    return message;
                  })
                  .map((error, index) => (
                    <Fragment key={index}>{error}</Fragment>
                  ))}
              </span>
            ) : null;

            /** field control to show */
            let control = <span></span>;

            /** ref function */
            const ref = (el) => {
              if (!el || !("setCustomValidity" in el)) return;
              /** also set native browser form validity to get various built-in features */
              if (isError) el?.setCustomValidity("Error with this field");
              else el?.setCustomValidity("");
            };

            /** single select */
            if (_enum?.length) {
              const options = uniq(["", ..._enum.filter(Boolean)]).map(
                (value) => ({
                  value,
                  label: upperFirst(value),
                }),
              );
              control = (
                <Select
                  options={options}
                  value={value || ""}
                  onChange={onChange}
                />
              );
            } else if (types.includes("string"))
              /** text input */
              control = (
                <TextBox
                  pattern={pattern}
                  value={value || ""}
                  onChange={(value) => onChange(value || null)}
                />
              );
            else if (types.includes("number") || types.includes("integer")) {
              /** min value allowed to be input */
              const min = Math.max(
                ...[minimum, data[greater_than]].filter(
                  (value) => typeof value === "number",
                ),
              );
              /** max value allowed to be input */
              const max = Math.min(
                ...[maximum, data[less_than]].filter(
                  (value) => typeof value === "number",
                ),
              );

              control = (
                <NumberBox
                  name={name}
                  label={label}
                  required={required}
                  pattern={pattern}
                  step={types.includes("integer") ? 1 : "any"}
                  min={Number.isFinite(min) ? min : undefined}
                  max={Number.isFinite(max) ? max : undefined}
                  value={value}
                  onChange={onChange}
                />
              );
            }

            return (
              <Fragment key={index}>
                {/* props that should be on all fields */}
                {cloneElement(control, { ref, name, label, required })}
                {error}
              </Fragment>
            );
          })}
        </section>
      ))}

      <section>
        <Button type="submit" design="bubble">
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </Form>
  );
};

export default SchemaForm;

const ObjectFieldTemplate = ({
  idSchema,
  schema,
  title,
  description,
  properties,
  formContext: { sections },
}) =>
  /** if at top level of schema */
  idSchema.$id === "root" ? (
    /** split fields into sections */
    sections.map((section) => (
      <section key={section}>
        {/* section top */}
        <Heading level={2} className={classes.full}>
          {section === undefined ? "Other" : section}
        </Heading>

        {/* section fields */}
        <div className={clsx(classes.grid, classes.card)}>
          {properties
            .filter(
              (property) =>
                schema.properties[property.name]?.section === section,
            )
            .map((property) => (
              <Fragment key={property.name}>{property.content}</Fragment>
            ))}
        </div>
      </section>
    ))
  ) : (
    /** otherwise, render children */
    <>
      <span className={classes.full}>{title}</span>
      <span className={classes.full}>{description}</span>
      {properties.map((property) => (
        <Fragment key={property.name}>{property.content}</Fragment>
      ))}
    </>
  );

const ArrayFieldTemplate = ({
  schema,
  items,
  formContext,
  canAdd,
  onAddClick,
}) => {
  const tooltip = getTooltip(schema);
  const label = getLabel(schema, tooltip);

  /** get different fields on property */
  const { enum: _enum } = schema;

  return (
    <div className={clsx("col", classes.card, classes.full)}>
      <span className={classes.heading}>{label}</span>

      <div className={classes.array}>
        {items.map(
          (
            {
              schema,
              children,
              onReorderClick,
              onDropIndexClick,
              onAddIndexClick,
            },
            index,
          ) => (
            <Fragment key={index}>
              <div className={clsx(classes.grid, classes.card)}>{children}</div>

              <div
                className={classes.actions}
                style={{
                  flexDirection:
                    schema.type === "object" &&
                    Object.keys(schema.properties)?.length >= 3
                      ? "column"
                      : "row",
                }}
              >
                <span>{index + 1}</span>
                <Button
                  data-tooltip="Move up"
                  style={{ gridColumn: "3" }}
                  onClick={formContext.triggerValidate(
                    onReorderClick(
                      index,
                      index > 0 ? index - 1 : items.length - 1,
                    ),
                  )}
                >
                  <FaArrowUp />
                </Button>
                <Button
                  data-tooltip="Move down"
                  onClick={formContext.triggerValidate(
                    onReorderClick(
                      index,
                      index < items.length - 1 ? index + 1 : 0,
                    ),
                  )}
                >
                  <FaArrowDown />
                </Button>
                <Button
                  data-tooltip="Remove"
                  onClick={formContext.triggerValidate(onDropIndexClick(index))}
                >
                  <FaTrash />
                </Button>
                <Button
                  data-tooltip="Insert new above"
                  onClick={formContext.triggerValidate(onAddIndexClick(index))}
                >
                  <FaPlus />
                </Button>
              </div>
            </Fragment>
          ),
        )}
      </div>

      <span className={classes.full}>
        {canAdd && (
          <Button
            data-tooltip="Add new"
            design="plain"
            onClick={formContext.triggerValidate(onAddClick)}
          >
            <FaPlus />
          </Button>
        )}
      </span>
    </div>
  );
};

/** add custom keyword for validating one value less than other value */
const lessThan = {
  id: "lessThan",
  keyword: "less_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    const fullData = node.context.data;
    const otherKey = node.schema.less_than;
    const otherValue = get(fullData, otherKey);
    if (thisValue <= otherValue) return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be less than ${otherKey}` },
    });
  },
};

/** add custom keyword for validating one value greater than other value */
const greaterThan = {
  id: "greaterThan",
  keyword: "greater_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    const fullData = node.context.data;
    const otherKey = node.schema.greater_than;
    const otherValue = get(fullData, otherKey);
    if (thisValue >= otherValue) return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be greater than ${otherKey}` },
    });
  },
};

/** https://www.npmjs.com/package/json-schema-library#user-content-extending-a-draft */
const draft = extendDraft(draft2020, {
  /** match all schema */
  $schemaRegEx: ".*",
  /** add custom validation rules */
  keywords: [lessThan, greaterThan],
});
