import { Fragment, useRef, useState } from "react";
import { FaArrowDown, FaArrowUp, FaPlus, FaTrash } from "react-icons/fa6";
import { LuSend } from "react-icons/lu";
import clsx from "clsx";
import { uniq, upperFirst } from "lodash-es";
import { useDebounceFn } from "@reactuses/core";
import Form from "@rjsf/core";
import Button from "@/components/Button";
import Heading from "@/components/Heading";
import Help from "@/components/Help";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TextBox from "@/components/TextBox";
import { validator } from "@/data/schema";
import classes from "./SchemaForm.module.css";

/** form automatically generated from json schema */
const SchemaForm = ({ schema, sections, data: initialData = {} }) => {
  const ref = useRef(null);
  const [formData, setFormData] = useState(initialData);

  /** revalidate whole form */
  const validate = useDebounceFn(() => ref.current?.validateForm(), 300);

  /** create func with debounced revalidation of whole form */
  /** https://github.com/rjsf-team/react-jsonschema-form/issues/3616#issuecomment-3221419763 */
  const triggerValidate =
    (func) =>
    (...args) => {
      func(...args);
      validate.run();
    };

  return (
    <Form
      ref={ref}
      className={classes.form}
      /** full form data */
      formData={formData}
      onChange={(event) => setFormData(event.formData)}
      /** allow certain info to be accessed from any field */
      formContext={{
        fullSchema: schema,
        fullData: formData,
        sections,
        triggerValidate,
      }}
      /** shape of data */
      schema={schema}
      /** custom validator */
      validator={validator}
      /** hide summary list of errors */
      showErrorList={false}
      /** never validate automatically, only on demand by calling validateForm */
      noValidate
      /** custom layouts */
      templates={{
        ObjectFieldTemplate,
        ArrayFieldTemplate,
        FieldTemplate,
        FieldErrorTemplate,
      }}
      /** custom fields */
      fields={{ StringField, NumberField }}
      onError={console.error}
      onSubmit={({ formData }) => console.info(formData)}
    >
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

/** @param {import("@rjsf/utils").ObjectFieldTemplateProps} props */
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

/** @param {import("@rjsf/utils").ArrayFieldTemplateProps} props */
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

/** @param {import("@rjsf/utils").FieldTemplateProps} props */
const FieldTemplate = ({ children, errors }) => (
  <>
    {children}
    {errors}
  </>
);

/** @param {import("@rjsf/utils").FieldErrorProps} props */
const FieldErrorTemplate = ({ errors }) =>
  errors?.length ? (
    <div className={classes.error}>
      {errors
        .map((error) => upperFirst(error))
        .join(", ")
        .replaceAll(/(\w),(\w)/g, "$1 | $2")}
    </div>
  ) : null;

/** @param {import("@rjsf/utils").WidgetProps} props */
const StringField = ({
  schema,
  name,
  formData: value,
  onChange,
  formContext: { fullSchema, triggerValidate },
}) => {
  const types = getTypes(schema);
  const required = getRequired(schema, fullSchema);
  const tooltip = getTooltip(schema);
  const label = getLabel(schema, tooltip);

  /** modify onChange func */
  onChange = triggerValidate(onChange);

  /** get different fields on property */
  const { pattern, enum: _enum } = schema;

  /** single select */
  if (_enum)
    return (
      <Select
        name={name}
        label={label}
        required={required}
        options={uniq(["", ..._enum.filter(Boolean)]).map((value) => ({
          value,
          label: upperFirst(value),
        }))}
        value={value || ""}
        onChange={onChange}
      />
    );

  /** text input */
  if (types.includes("string"))
    return (
      <TextBox
        name={name}
        label={label}
        required={required}
        pattern={pattern}
        value={value || ""}
        onChange={onChange}
      />
    );
};

/** @param {import("@rjsf/utils").WidgetProps} props */
const NumberField = ({
  schema,
  name,
  formData: value,
  onChange,
  formContext: { fullSchema, fullData, triggerValidate },
}) => {
  const types = getTypes(schema);
  const required = getRequired(schema, fullSchema);
  const tooltip = getTooltip(schema);
  const label = getLabel(schema, tooltip);

  /** modify onChange func */
  onChange = triggerValidate(onChange);

  /** get different fields on property */
  const {
    pattern,
    enum: _enum,
    minimum,
    maximum,
    less_than,
    greater_than,
  } = schema;

  /** min value allowed to be input */
  const min = Math.max(
    ...[minimum, fullData[greater_than]].filter(
      (value) => typeof value === "number",
    ),
  );
  /** max value allowed to be input */
  const max = Math.min(
    ...[maximum, fullData[less_than]].filter(
      (value) => typeof value === "number",
    ),
  );

  return (
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
};

/** normalize type to array */
const getTypes = ({ type }) => [type].flat();

/** make help tooltip to show */
const getTooltip = ({ description, examples }) =>
  [
    description,
    ...(examples?.length ? [" "] : []),
    examples?.length ? `Examples: ${examples.join(", ")}` : "",
  ]
    .filter(Boolean)
    .join("<br/>");

/** is the field required to be set */
const getRequired = ({ name, type }, { required }) =>
  required.includes(name) && !getTypes(type).includes("null");

/** full label elements to show */
const getLabel = ({ title, name }, tooltip) => (
  <>
    {title ?? name ?? ""}
    {tooltip && <Help>{tooltip}</Help>}
  </>
);
