import { cloneElement, Fragment, useMemo, useState } from "react";
import { FaArrowDown, FaArrowUp, FaPlus, FaTrash } from "react-icons/fa6";
import { LuSend } from "react-icons/lu";
import clsx from "clsx";
import { compileSchema, draft2020, extendDraft } from "json-schema-library";
import { cloneDeep, uniq, upperFirst } from "lodash-es";
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
const SchemaForm = ({ schema, sections, data: initialData }) => {
  /** compile schema */
  const rootNode = useMemo(
    () => compileSchema(schema, { drafts: [draft] }),
    [schema],
  );

  /** if no initial data, generate blank values */
  initialData ??= rootNode.getData(null, { extendDefaults: false });

  /** current form data state */
  const [data, setData] = useState(initialData);

  /** revalidate when data changes */
  const { errors } = useMemo(() => {
    rootNode.context.data = data;
    return rootNode.validate(data);
  }, [rootNode, data]);

  return (
    <Form className={classes.form} onSubmit={console.info}>
      {/* correct section coloring */}
      <span></span>

      {sections.map((section, index) => (
        <section key={index}>
          <Heading level={2}>{section}</Heading>

          <Field
            rootNode={rootNode}
            schema={schema}
            section={section}
            node={rootNode}
            data={data}
            setData={setData}
            errors={errors}
          />
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

/**
 * @param {Object} props
 * @param {import("json-schema-library").SchemaNode} props.node
 * @param {import("json-schema-library").JsonError[]} props.errors
 */
const Field = ({
  rootNode,
  schema,
  section,
  index,
  node,
  data,
  setData,
  errors,
}) => {
  /** path to field, without root part */
  let path = node.evaluationPath.replace(/^#\/properties/, "");

  if (index !== undefined) path = path.replace(/\/items$/, `/${index}`);

  /** are we at top level of schema */
  const level = split(path).length;

  /** filter out fields that should not be displayed in this section */
  if (level === 1 && node.schema.section !== section) return;

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
  } = node.schema;

  /** normalize type to array */
  const types = [type].flat();

  /** form data name */
  const name = split(path).pop();

  /** help tooltip to show */
  const tooltip = [
    description,
    ...(examples?.length ? [" "] : []),
    examples?.length ? `Examples: ${examples.join(", ")}` : "",
  ]
    .filter(Boolean)
    .join("<br/>");

  /** is the field required to be set */
  const required = schema.required.includes(name) && !types.includes("null");

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
  const onChange = (value) => {
    console.log(path);
    setData(set(cloneDeep(data), path, value));
  };

  /** get validation errors associated with this field */
  const fieldErrors = errors.filter(({ data }) => data.pointer.endsWith(path));
  /** is there an error on this field */
  const isError = !!fieldErrors.length;

  /** error component to show */
  const error = isError ? (
    <span className={clsx(classes.error, classes.full)}>
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
          if (code === "compare-value-error") return <>{data.schema.message}</>;
          if (code === "type-error")
            return <>Must be {makeList(data.expected, "code")}</>;
          return message;
        })
        .map((error, index) => (
          <Fragment key={index}>{error}</Fragment>
        ))}
    </span>
  ) : null;

  /** field control to show */
  let control = <span className={classes.full}>missing control</span>;

  /** ref function */
  const ref = (el) => {
    if (!el || !("setCustomValidity" in el)) return;
    /** also set native browser form validity to get various built-in features */
    if (isError) el?.setCustomValidity("Error with this field");
    else el?.setCustomValidity("");
  };

  if (types.includes("object"))
    /** object group */
    control = (
      <div className={clsx(classes.full, classes.grid, classes.card)}>
        {level > 0 && (
          <div className={clsx(classes.full, classes.heading)}>{label}</div>
        )}
        {Object.keys(node.schema.properties).map((key) => {
          return (
            <Field
              key={key}
              rootNode={rootNode}
              schema={schema}
              section={section}
              node={node.getChildSelection(key)[0]}
              data={data}
              setData={setData}
              errors={errors}
            />
          );
        })}
      </div>
    );
  else if (types.includes("array")) {
    /** array group */
    const items = get(data, path) ?? [];
    control = (
      <div className={clsx(classes.full, classes.grid, classes.card)}>
        <div className={clsx(classes.full, classes.heading)}>{label}</div>
        <div className={clsx(classes.full, classes.array)}>
          {items.map((item, index) => {
            return (
              <Fragment key={index}>
                <Field
                  key={index}
                  index={index}
                  schema={schema}
                  section={section}
                  node={node.getChildSelection("items")[0]}
                  data={data}
                  setData={setData}
                  errors={errors}
                />
                <div className={classes.actions}>
                  <span>{index + 1}</span>
                  <Button data-tooltip="Move up">
                    <FaArrowUp />
                  </Button>
                  <Button data-tooltip="Move down">
                    <FaArrowDown />
                  </Button>
                  <Button data-tooltip="Remove">
                    <FaTrash />
                  </Button>
                </div>
              </Fragment>
            );
          })}
        </div>
        <Button
          design="plain"
          onClick={() => {
            const _data = cloneDeep(data);
            const list = get(_data, path) ?? [];
            let newItem = node
              .getNodeChild("items")
              .node.getData(undefined, {});
            if (newItem === "") newItem = null;
            setData(set(_data, path, [...list, newItem]));
          }}
        >
          <FaPlus />
          <span>Add</span>
        </Button>
      </div>
    );
  } else if (_enum?.length) {
    /** single select */
    const options = uniq(["", ..._enum.filter(Boolean)]).map((value) => ({
      value,
      label: upperFirst(value),
    }));
    control = (
      <Select options={options} value={value || ""} onChange={onChange} />
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
    <Fragment>
      {/* props that should be on all fields */}
      {cloneElement(control, { ref, name, label, required })}
      {error}
    </Fragment>
  );
};
